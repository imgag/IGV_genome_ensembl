"""
    Generates a IGV genome file in the new JSON format from a template (including updating gene file)
"""
import argparse
import json
import operator
import os
import subprocess
import tempfile
import urllib.request
import gzip
import shutil

# global variables
genome_json = {}

"""
sorting order for non integer chromosomes
"""
sorting_order = {'X': 100, 'Y': 101, 'MT': 102}

def parse_args():
    """
                parses the arguments

    :return:    argparse object containing all provided arguments
    """

    print("parsing args...")

    parser = argparse.ArgumentParser(description="Generates a IGV genome JSON")
    parser.add_argument("template_file", help="Template JSON containing all links to the input files.")
    parser.add_argument("hgnc_file", help="file path to the the HGNC table file (containing HGNC id <-> gene name mapping")
    parser.add_argument("output", help="file path for the generated output IGV genome JSON file")

    return parser.parse_args()


def parse_json(template_file: str):
    """
            parses the template file

    :param template_file:   file path to the json template file
    """
    global genome_json
    with open(template_file, 'r') as json_file:
        genome_json = json.load(json_file)
    return


def download_files(output_folder: str):
    """
            downloads all distant files in the template json and links to them

    :param output_folder:   target folder for the downloads

    :return:                modified JSON template with links to the local files
    """
    global genome_json
    # download all required files

    for key in ["fastaURL", "indexURL", "cytobandURL", "aliasURL"]:
        url = genome_json[key]
        filename = os.path.basename(url)
        print("Downloading file '" + filename + "'...")
        urllib.request.urlretrieve(url, os.path.join(output_folder, filename))
        genome_json[key] = filename

    # download tracks
    for track in genome_json["tracks"]:
        url = track["url"]
        filename = os.path.basename(url)
        print("Downloading file '" + filename + "'...")
        urllib.request.urlretrieve(url, os.path.join(output_folder, filename))
        track["url"] = filename
        # load optional index
        if "indexURL" in track:
            url = track["indexURL"]
            filename = os.path.basename(url)
            print("Downloading file '" + filename + "'...")
            urllib.request.urlretrieve(url, os.path.join(output_folder, filename))
            track["indexURL"] = filename

    return


def load_hgnc_file(hgnc_filepath):
    print("Parsing HGNC file...")
    hgnc_mapping = {}
    with open(hgnc_filepath, 'r', encoding="utf8") as hgnc_file:
        for line in hgnc_file:
            if line.startswith("hgnc_id\tsymbol\tname"):
                # skip header
                continue
            elif line.startswith("HGNC:"):
                # parse line
                split_line = line.split('\t')
                hgnc_id = int(split_line[0].split(':')[1])
                symbol = split_line[1].strip()
                hgnc_mapping[hgnc_id] = symbol
    print("\t " + str(len(hgnc_mapping.keys())) + " ids parsed.")
    return hgnc_mapping


def sort_gff3(content):
    """
                sorts the content gff3 file based on the sorting order defined above

    :param content:     content of the gff file (without headers)
    :return:            sorted list of lists
    """

    print("sorting gff3 data...")

    # preprocess chromosome column and positioning column
    for line in content:
        if line[0] in sorting_order:
            line[0] = sorting_order[line[0]]
        else:
            line[0] = int(line[0])

        line[3] = int(line[3])
        line[4] = int(line[4])

    # sort content
    content.sort(key=operator.itemgetter(0, 3))

    # create reverse sorting_order:
    reverse_sorting_order = {}
    for k, v in sorting_order.items():
        reverse_sorting_order[v] = k

    # replace placeholder with original chromosome names:
    for line in content:
        if line[0] in reverse_sorting_order:
            line[0] = reverse_sorting_order[line[0]]

    return content

def update_gene_file(hgnc_mapping, output_folder):
    """

    :param hgnc_mapping:
    :return:
    """

    print("Modifying GFF3 files (updating gene names) ...")
    global genome_json
    # modify all gff3 track files:
    for track in genome_json["tracks"]:
        if "format" in track and track["format"] == "gff3":

            print("Modifying GFF3 file '" + track["url"] + "'...")

            # unzip and modify
            n_comment_lines = 0
            n_unmodified_lines = 0
            n_modified_lines = 0
            n_ignored = 0
            comment_lines = []
            content_lines = []
            with gzip.open(os.path.join(output_folder, track["url"]), 'rb') as compressed_gff3:

                    for line in compressed_gff3:
                        line = line.decode('utf-8')
                        # skip comments
                        if line.startswith("#"):
                            if line.strip() == "###":
                                # ignore
                                n_ignored += 1
                                continue
                            comment_lines.append(line)
                            n_comment_lines += 1
                            continue
                        # detect HGNC ids
                        annotation_column = line.split('\t')[8]
                        if "[Source:HGNC Symbol%3BAcc:" in annotation_column:
                            kv_list = annotation_column.split(';')
                            idx_name = -1
                            idx_description = -1
                            for idx in range(len(kv_list)):
                                if kv_list[idx].startswith("Name="):
                                    idx_name = idx
                                elif kv_list[idx].startswith("description="):
                                    idx_description = idx
                            if idx_description > -1:
                                # extract HGNC id
                                hgnc_id = int(kv_list[idx_description].split('[')[1].split(']')[0].split(':')[-1])

                                if hgnc_id not in hgnc_mapping.keys():
                                    print("Warning: HGNC id " + str(hgnc_id) + " not found in HGNC file!")
                                    # store line unmodified
                                    content_lines.append(line.split('\t'))
                                    n_unmodified_lines += 1
                                    continue
                                if idx_name > -1:
                                    kv_list[idx_name] = "Name=" + hgnc_mapping[hgnc_id]
                                else:
                                    kv_list.append("Name=" + hgnc_mapping[hgnc_id])

                                split_line = line.split('\t')
                                split_line[8] = ";".join(kv_list)
                                content_lines.append(split_line)
                                n_modified_lines += 1
                        else:
                            # no HGNC identifier -> store unmodified
                            content_lines.append(line.split('\t'))
                            n_unmodified_lines += 1

            # stats
            print("\tcomment lines: " + str(n_comment_lines))
            print("\tunmodified lines: " + str(n_unmodified_lines))
            print("\tmodified lines: " + str(n_modified_lines))
            print("\tignored lines: " + str(n_ignored))

            # sort content

            content_lines = sort_gff3(content_lines)

            # write modified file to disk
            print("Writing modified file to disk...")
            with open(os.path.join(output_folder, os.path.splitext(track["url"])[0]), 'wt') as modified_gff3:
                for line in comment_lines:
                    modified_gff3.write(line)
                for line in content_lines:
                    modified_gff3.write("\t".join(str(e) for e in line))

            # bgzip
            print("Compressing file...")
            print("\tignored lines: " + str(n_ignored))
            rc = subprocess.call(["bgzip", "-f", os.path.join(output_folder, os.path.splitext(track["url"])[0])])
            if rc != 0:
                raise RuntimeError("bgzip failed with return code " + str(rc) + "!")

            # tabix
            print("Indexing file...")
            rc = subprocess.call(["tabix", "-p", "gff", os.path.join(output_folder, track["url"])])
            if rc != 0:
                raise RuntimeError("tabix failed with return code " + str(rc) + "!")
            # add index to JSON
            track["indexURL"] = track["url"] + ".tbi"
    return


def update_alias_file(output_folder):
    """
                Extends the alias tab file
    """
    print("Updating alias file...")
    alias_file_name = genome_json["aliasURL"]
    file_buffer = []
    with open(os.path.join(output_folder, alias_file_name), 'r') as alias_file:
        for line in alias_file:
            if line.startswith("chrM"):
                # add 'chrMT' and 'M' as valid aliases
                line = line.strip() + "\tchrMT\tM\n"
            file_buffer.append(line)

    with open(os.path.join(output_folder, alias_file_name), 'w') as alias_file:
        alias_file.writelines(file_buffer)
    return


def main():

    args = parse_args()

    # read template
    parse_json(args.template_file)

    # download files to local storage
    output_folder = os.path.dirname(args.output)
    download_files(output_folder)

    # load hgnc file
    hgnc_mapping = load_hgnc_file(args.hgnc_file)

    # update gene file
    update_gene_file(hgnc_mapping, output_folder)

    # update alias file
    update_alias_file(output_folder)

    # store modified JSON file
    with open(args.output, 'w') as output_file:
        print("Writing genome JSON file...")
        json.dump(genome_json, output_file, indent=4)

    print("\nfinished.")


if __name__ == '__main__':
    main()
