package dvec

// Fragment is a single building block stretch of DNA for assembly
type Fragment struct {
	// ID is a unique identifier for this fragment
	ID string

	// fragment's sequence (linear)
	Seq string

	// primers necessary to create this (if pcr fragment)
	Primers []Primer

	// Entry of this fragment In the DB that it came from
	// Used to look for off-targets
	Entry string
}

// Match is a blast "hit" in the blastdb
type Match struct {
	// ID of the matched fragment in the database
	ID string

	// Seq of the match on the target vector
	Seq string

	// Start of the fragment (0-indexed)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Whether it's a circular fragment (vector, plasmid, etc)
	Circular bool

	// The number of mismatching bps in the match (for primer off-targets)
	Mismatch int
}

// Length returns the length of the match on the target fragment
func (m *Match) Length() int {
	return m.End - m.Start + 1 // it's inclusive
}

// Primer is a single Primer used to ampligy a parent fragment
type Primer struct {
	// Seq of the primer (In 5' to 3' direction)
	Seq string

	// Strand of the primer; true if template, false if complement
	Strand bool

	// Start of the fragment (0-index)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Penalty score
	Penalty float32

	// PairPenalty score from primer3
	PairPenalty float32

	// Tm of the primer
	Tm float32

	// GC % max
	GC float32
}