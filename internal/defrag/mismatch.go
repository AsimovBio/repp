package defrag

import (
	"fmt"
	"io/ioutil"
	"os/exec"
	"path"
	"path/filepath"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// isMismatch reutns whether the match constitutes a mismatch
// between it and the would be primer sequence
//
// source: http://depts.washington.edu/bakerpg/primertemp/
//
// The equation used for the melting temperature is:
// Tm = 81.5 + 0.41(%GC) - 675/N - % mismatch, where N = total number of bases.
func isMismatch(match Match) bool {
	primer := strings.ToLower(match.Seq)
	primerL := float32(len(primer))

	noA := strings.Replace(primer, "a", "", -1)
	noT := strings.Replace(noA, "t", "", -1)
	gcPerc := float32(len(noT)) / primerL
	tmNoMismatch := 81.5 + 0.41*gcPerc - 675/float32(len(primer))
	tmWithMismatch := tmNoMismatch - float32(match.Mismatch)/primerL

	return tmWithMismatch > 40 // TODO: move to settings
}

// Mismatch finds mismatching sequences between the query sequence and
// the parent sequence
//
// The parent sequence is passed as the entry id as it exists in the blast db
// db is passed as the path to the db we're blasting against
func Mismatch(primer, parent, db string, v config.VendorConfig) (mismatch bool, match Match, err error) {
	// path to the entry batch file to hold the parent entry accession
	entry, _ := filepath.Abs(path.Join(v.Blastdir, parent+".entry"))

	// path to the output sequence file  from querying the parent's sequence from the BLAST db
	parentPath, _ := filepath.Abs(path.Join(v.Blastdir, parent+".out"))

	// path the query sequence input file
	in, _ := filepath.Abs(path.Join(v.Blastdir, parent+".primer.query"))

	// path to the blastOutput file
	out, _ := filepath.Abs(path.Join(v.Blastdir, parent+".blast"))

	// write entry to file
	// this was a 2-day bug I couldn't resolve...
	// I was using the "-entry" flag on exec.Command, but have since
	// switched to the simpler -entry_batch command (on a file) that resolves the issue
	if err = ioutil.WriteFile(entry, []byte(parent), 0666); err != nil {
		return false, match, fmt.Errorf("Failed to write batch entry list: %v", err)
	}

	// make a blastdbcmd command (for querying a DB, very different from blastn)
	queryCmd := exec.Command(
		v.Makeblastdb,
		"-db", db,
		"-dbtype", "nucl",
		"-entry_batch", entry,
		"-out", parentPath,
		"-outfmt", "%f", // fasta format
	)
	if _, err := queryCmd.CombinedOutput(); err != nil {
		return false, match, nil // pretending there wasn't any mismatch
	}

	// create blast input file
	inContent := fmt.Sprintf(">primer\n%s\n", primer)
	if err = ioutil.WriteFile(in, []byte(inContent), 0666); err != nil {
		return false, match, fmt.Errorf("Failed to write primer sequence to query FASTA file: %v", err)
	}

	// blast the query sequence against the parent sequence
	b := blastExec{
		in:      in,
		out:     out,
		subject: parentPath,
		blastn:  v.Blastn,
	}

	// execute blast
	if err = b.runAgainst(); err != nil {
		return false, match, fmt.Errorf("Failed to run blast against parent: %v", err)
	}

	// get the BLAST matches
	matches, err := b.parse()
	if err != nil {
		return false, Match{}, fmt.Errorf("Failed to parse matches from %s: %v", out, err)
	}

	// parse the results and check whether any are cause for concern (by Tm)
	primerCount := 1 // times we expect to see the primer itself
	for i, m := range matches {
		if i == 0 && m.Circular {
			// if the match is against a circular fragment, we might expect to see
			// the primer's sequence twice, rather than just once
			primerCount++
		}

		// one of the matches will, of course, be against the primer itself
		// and we don't want to double count it
		if primerCount > 0 && m.Seq == primer {
			primerCount--
			continue
		} else if isMismatch(m) {
			return true, m, nil
		}
	}

	return false, Match{}, nil
}