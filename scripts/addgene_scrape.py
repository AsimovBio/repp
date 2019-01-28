
import multiprocessing
import os
import string
import sys
from concurrent.futures import ThreadPoolExecutor, wait

import requests
from lxml import html

# output test directory for addgene
os.chdir(os.path.join("assets", "addgene"))

# valid file characters
VALID_CHARS = "-_.() %s%s" % (string.ascii_letters, string.digits)


def clean_name(addgene_name):
    """convert an addgene name to one that's file name compliant
    """
    return "".join([c for c in addgene_name if c in VALID_CHARS])


"""
Scrape the AddGene database into the local test/addgene directory for development/testing

AddGene pages are all indexed by the vector's id like: https://www.addgene.org/8389/sequences/

create FASTA files for each in test/data/addgene

if Addgene has full sequence information (prefered):
    save it as ">gnl|addgene|<addgene id> name-circular" in the fasta file
else:
    save it as ">gnl|addgene|<addgene id>.<linear index> name-linear" in the fasta file

save the file as <vector_id>_<addgene id>.fa
"""


def parse(page_index):
    """get a page from addgene and save its fragments/vectors if it's a defined id
    """
    title = '//*[@id="page-body"]/main/div/h1/span/text()'
    sequence_text_area = '//*[@id="addgene-full"]/ul/li/div/div/div[2]/textarea/text()'
    partial_sequence = '//*[@id="addgene-partial"]/ul/li/div/div/div[2]/textarea/text()'

    page = requests.get("https://www.addgene.org/" + str(page_index))
    if "Discontinued" in str(page.content):
        # avoid discontinued sequences
        return False
    page = requests.get("https://www.addgene.org/" + str(page_index) + "/sequences/")
    tree = html.fromstring(page.content)

    # XPath to the title

    name = tree.xpath(title)
    if not name:  # many ids are unused/not on the site
        return
    name = clean_name(name[0])

    # XPath to a seq block for full sequence information
    seq_blocks = tree.xpath(sequence_text_area)

    # create the file, we'll delete if afterwards if there was only
    # partial DNA available (had N's in it)
    file_path = os.path.join("files", str(page_index) + "_" + name + ".fa")
    file_written = False
    with open(file_path, "w") as f:
        index = str(page_index)

        if seq_blocks:
            # full seq information available, write that
            seq = seq_blocks[0].split("\n")[1:]
            seq = "".join("".join(seq).split())

            out = ">gnl|addgene|{} {}-circular\n{}\n".format(index, name, seq + seq)
            f.write(out)
            file_written = True
        else:
            # partial seq information available, write that
            partial_seq_blocks = tree.xpath(partial_sequence)

            count = 1
            for partial_seq in partial_seq_blocks:
                seq = partial_seq.split("\n")[1:]
                seq = "".join("".join(seq).split())
                frag_count = str(count)

                # don't save files with unknown basepairs
                if "N" in seq:
                    continue

                out = ">gnl|addgene|{}.{} {}-linear\n{}\n".format(
                    index, frag_count, name, seq
                )
                count += 1
                file_written = True

    if not file_written:  # was all partial data, remove
        os.remove(file_path)

    return file_written


assert parse(39081) == True
assert parse(24873) == False  # discontinued


def scrape_files():
    INDEX_LIMIT = 150000

    # find the saved files and Addgene pages already seen on local filesystem
    files = [f for f in os.listdir(".") if f and "fa" in f]
    max_page = max([int(f.split("_")[0]) for f in files] + [1])

    # files to ignore, don't correspond to addgene parts
    addgeneignore = set(
        [int(f.strip()) for f in open(".addgeneignore", "r").readlines()]
    )

    page_indicies = [
        i for i in range(1, INDEX_LIMIT) if i > max_page and i not in addgeneignore
    ]

    futures = []
    with ThreadPoolExecutor() as executor:
        """try to scrape with all available CPUs (it's x5 off that since it's IO dependent)
        """
        for page in page_indicies:
            futures.append(executor.submit(parse, page))
    wait(futures)

    # recalculate the seen indicies and save to the fs so we don't scan again in future
    seen_indices = set(
        [int(f.split("_")[0]) for f in os.listdir(".") if "_" in f and f is not ".DS"]
    )
    empty_indices = [i for i in range(1, INDEX_LIMIT) if i not in seen_indices]
    with open(".addgeneignore", "w") as f:
        for ind in empty_indices:
            f.write(str(ind) + "\n")


def combine():
    # combine all the FASTA files into one
    os.chdir("db")

    # open a single FASTA for a combined db
    with open("addgene", "w") as combined_fasta:
        os.chdir(os.path.join("..", "files"))

        for f in os.listdir("."):
            file_contents = open(f, "r").read().strip() + "\n"
            combined_fasta.write(file_contents)


# scrape_files()
# combine()
