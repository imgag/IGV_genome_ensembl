"""
Converts a gff3 file to genePred
"""
import operator
import shutil
import stat
import subprocess
import sys
import argparse
import os
import tempfile
from shutil import copyfile
from six.moves.urllib.request import urlopen
from six.moves.urllib.error import URLError, HTTPError
import zipfile
from collections import OrderedDict

"""
sorting order for non integer chromosomes
"""
sorting_order = {'X': 100, 'Y': 101, 'MT': 102}

"""
download URL to the gff3ToGenePred tool 
"""
gff_converter_url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred"


def parse_args():
    """
                parses the arguments

    :return:    argparse object containing all provided arguments
    """

    print("parsing args...")

    parser = argparse.ArgumentParser(description="Converts gff3 files into genePred format")
    parser.add_argument("gff_file", help="file path to the gff3 input file (with ensembl annotation)")
    parser.add_argument("hgnc_file",
                        help="file path to the the HGNC table file (containing HGNC id <-> gene name mapping")
    parser.add_argument("genome_file", help="file path to a IGV .genome file which is then used to generate own .genome"
                                            + " file with the generated annotation file")
    parser.add_argument("output", help="file path for the generated output IGV .genome file")

    return parser.parse_args()


def validate_args(args):
    """
                validates all given parameters

    :param args:    argparse object containing all given arguments
    :return:        True if parameter are correct, else False
    """

    print("validating args...")

    # check gff file
    if not os.path.isfile(args.gff_file):
        sys.stderr.write("gff file not found in %s" % args.gff)
        return False
    # check HGNC file
    if not os.path.isfile(args.hgnc_file):
        sys.stderr.write("hgnc mapping file not found in %s" % args.gff)
        return False
    # check genome file
    if not os.path.isfile(args.genome_file):
        sys.stderr.write("genome file not found in %s" % args.gff)
        return False

    return True


def read_gff_file(gff_file_path):
    """
                reads the gff file and returns it divided in header and content

    :param gff_file_path:   file path to the gff file
    :return:
        header              list of all header lines
        content             content of the gff file as list of lists
    """

    print("reading gff file...")

    with open(gff_file_path, 'r') as gff_file:
        data = gff_file.readlines()

    header = []
    content = []

    replaced_genes = 0
    replaced_transcripts = 0

    for line in data:
        if line.startswith('#'):
            if not line.startswith('###'):
                header.append(line.strip())
        else:
            # ignore GL000xxx entries:
            if line.startswith('GL000'):
                continue
            # ignore KI270xxx entries:
            if line.startswith('KI270'):
                continue

            # split line by tab
            split_line = line.strip().split('\t')

            # treat all entries with ENST id as transcript and entries with ENSG id as genes
            # generate temp dict with all values:
            meta_data_dict = {entry.split('=')[0]: entry.split('=')[1] for entry in split_line[8].split(';')}
            if "ID" in meta_data_dict.keys():
                if meta_data_dict["ID"].startswith("transcript:ENST") and split_line[2] != "transcript":
                    split_line[2] = "transcript"
                    replaced_transcripts += 1
                elif meta_data_dict["ID"].startswith("gene:ENSG") and split_line[2] != "gene" and \
                        not split_line[2].endswith("_gene_segment"):
                    split_line[2] = "gene"
                    replaced_genes += 1

            content.append(split_line)

    # report replaced genes and transcripts
    print("\t {} entries with ENSG ids are treated as genes".format(replaced_genes))
    print("\t {} entries with ENST ids are treated as transcripts".format(replaced_transcripts))

    return header, content


def sort_gff(content):
    """
                sorts the content gff file based on the sorting order defined above

    :param content:     content of the gff file (without headers)
    :return:            sorted list of lists
    """

    print("sorting gff data...")

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


def write_gff(gff_file_handle, header, content):
    """
                writes the header and content to a gff file

    :param gff_file_handle:     file handle for the gff file
    :param header:              header lines (starting with #)
    :param content:             content of the gff file (as list of lists)
    :return:
    """

    print("writing gff file...")

    # write header:
    for line in header:
        gff_file_handle.write(line)
        gff_file_handle.write("\n")

    # write content:
    for line in content:
        gff_file_handle.write("\t".join(str(e) for e in line))
        gff_file_handle.write("\n")

    return


def setup_gff_converter(file_handle):
    """
                downloads the gff converter tool from the ucsc website

    :param file_handle:      file handle to save the downloaded binary
    :return:
    """

    print("initializing gff3ToGenePred converter...")

    # download gff converter
    try:
        print("downloading converter from: " + gff_converter_url)
        f = urlopen(gff_converter_url)

        # write downloaded file to disk
        file_handle.write(f.read())

    # handle errors
    except HTTPError as e:
        print("HTTP Error:", e.code, gff_converter_url)
    except URLError as e:
        print("URL Error:", e.reason, gff_converter_url)

    # make tool executable:
    os.chmod(file_handle.name, os.stat(file_handle.name).st_mode | stat.S_IEXEC)


def run_gff_converter(converter_file_handle, gff_file_handle, gene_pred_file_handle):
    """
                converts a gff3 file into the genePred format using the gff3ToGenePred tool from the ucsc website

    :param converter_file_handle:     file handle for the gff3ToGenePred binary
    :param gff_file_handle:           file handle for the input gff3 file
    :param gene_pred_file_handle:     file handle for the genePred output file

    :return: gene_pred_file_handle:   containing the modified and reopened tempfile
    """

    print("converting gff file...")

    # close temporary files before converting
    converter_file_handle.close()
    gff_file_handle.close()
    gene_pred_file_handle.close()

    # concat the command
    cmd = [converter_file_handle.name, gff_file_handle.name, gene_pred_file_handle.name]
    subprocess.call(cmd)

    # delete sorted gff file and converter (not required anymore)
    os.remove(gff_file_handle.name)
    os.remove(converter_file_handle.name)

    return open(gene_pred_file_handle.name, "r")


def read_hgnc_file(file_path):
    """
                reads a hgnc file and returns two dictionaries gene->id and id->gene

    :param file_path:        file path to the HGNC file
    :return:
            gene_to_hgnc:   dict mapping genes to HGNC ids
            hgnc_to_gene:   dict mapping HGNC ids to gene names
    """

    print("reading HGNC file...")

    # read hgnc file
    with open(file_path, 'r') as hgnc_file:
        hgnc_data = hgnc_file.readlines()

    # process data
    gene_to_hgnc = {}
    hgnc_to_gene = {}

    # start at 1 to skip header
    for line in hgnc_data[1:]:
        split_line = line.split('\t')
        hgnc_id = int(split_line[0].split(':')[1])
        gene_name = split_line[1]

        # store in dicts
        gene_to_hgnc[gene_name] = hgnc_id
        hgnc_to_gene[hgnc_id] = gene_name

    print("\t %i HGNC ids read" % len(hgnc_to_gene))

    return gene_to_hgnc, hgnc_to_gene


def generate_ensg_hgnc_mapping(gff_file_path, valid_hgnc_ids):
    """
                extracts a ENSG<->HGNC mapping from the gff3 file

    :param gff_file_path:   file path to the gff3 file
    :param valid_hgnc_ids:  set with valid HGNC ids
    :return:
            ensg_to_hgnc:   dict mapping ensembl gene id to HGNC ids
            hgnc_to_ensg:   dict mapping HGNC ids to ensembl gene id
    """

    print("generating ensembl gene id <-> HGNC mapping...")

    ensg_to_hgnc = {}
    hgnc_to_ensg = {}

    ensg_to_non_hgnc_gene = {}
    non_hgnc_gene_to_ensg = {}

    ignored_genes = []
    ignored_genes_dots = []

    # read file content without header
    with open(gff_file_path, 'r') as gff_file:
        gff_data = [line for line in gff_file.readlines() if not line.startswith('#')]

    # remove all non gene lines and split by tab
    gff_gene_lines = [line.split('\t') for line in gff_data if
                      line.split('\t')[2] in ["gene", "pseudogene", "processed_transcript", "RNA"]
                      or line.split('\t')[2].endswith('_gene_segment') or line.split('\t')[2].endswith('_gene')]

    # extract ENSG<->HGNC mapping:
    for line in gff_gene_lines:
        meta_data = line[8].split(';')
        # generate temp dict with all values:
        meta_data_dict = {entry.split('=')[0]: entry.split('=')[1] for entry in meta_data}

        # skip entries which do not have a ensembl gene id:
        try:
            ensg = meta_data_dict["gene_id"]
        except KeyError:
            continue

        # tries to extract HGNC id or saves gene name instead
        try:
            # skip all non HGNC ids:
            if "Source:HGNC Symbol%3BAcc:" not in meta_data_dict["description"]:
                raise ValueError
            hgnc = int(meta_data_dict["description"].split(':')[-1][:-1])
            # skip all invalid HGNC ids
            if hgnc not in valid_hgnc_ids:
                raise ValueError
        except (ValueError, KeyError):
            if '.' in meta_data_dict["Name"] or '-' in meta_data_dict["Name"]:
                ignored_genes_dots.append(meta_data_dict["Name"])
            else:
                ignored_genes.append(meta_data_dict["Name"])

            # use gene name instead of HGNC id:
            ensg_to_non_hgnc_gene[ensg] = meta_data_dict["Name"]
            non_hgnc_gene_to_ensg[meta_data_dict["Name"]] = ensg
            continue

        ensg_to_hgnc[ensg] = hgnc
        hgnc_to_ensg[hgnc] = ensg

    print("\t %i genes with HGNC ids mapped" % len(ensg_to_hgnc))
    print("\t %i genes without HGNC ids mapped" % len(ensg_to_non_hgnc_gene))

    return ensg_to_hgnc, hgnc_to_ensg, ensg_to_non_hgnc_gene, non_hgnc_gene_to_ensg


def read_gene_pred_file(gene_pred_file_handle):
    """
                reads a genePred file as list of lists

    :param gene_pred_file_handle:   file handle for the genePred input file
    :return:                        list of lists of all entries in the file
    """

    print("reading genePred file...")

    raw_data = gene_pred_file_handle.readlines()

    gene_pred_table = [line.split('\t') for line in raw_data]

    # close and delete file
    print("\t" + gene_pred_file_handle.name)
    gene_pred_file_handle.close()
    os.remove(gene_pred_file_handle.name)

    return gene_pred_table


def modify_gene_pred_data(gene_pred_data, ensg_to_hgnc, hgnc_to_gene, ensg_to_non_hgnc_gene):
    """
                modifies the genePred data to match the IGV requirements using HGNC ids

    :param gene_pred_data:          content of the genePred file as list of lists
    :param ensg_to_hgnc:            dict with mapping ensembl id -> HGNC id
    :param hgnc_to_gene:            dict with mapping HGNC id -> gene name
    :param ensg_to_non_hgnc_gene:   dict with mapping ensembl id -> gene name (which do not have a HGNC id)
    :return:                        list of lists with all entries of the modified genePred data
    """

    print("modifying genePred file...")

    hgnc_genes = 0
    non_hgnc_genes = 0
    for line in gene_pred_data:

        # remove "transcript:" in front of the ENST id:
        line[0] = line[0].split(':')[1]

        # get ENSG id
        ensg = line[11].split(':')[1].strip()

        # replace ENSG id with gene name
        try:
            gene_name = hgnc_to_gene[ensg_to_hgnc[ensg]]
            hgnc_genes += 1
        except KeyError:
            gene_name = ensg_to_non_hgnc_gene[ensg]
            non_hgnc_genes += 1
        line[11] = gene_name

        # add ENSG number as gene id in the first column
        ensg_int = int(ensg[4:])
        line.insert(0, str(ensg_int))

    print("\t gene names from %i hgnc genes and %i non hgnc genes were replaced" % (hgnc_genes, non_hgnc_genes))

    return gene_pred_data


def write_gene_pred_file(gene_pred_data, gene_pred_file_handle):
    """
                writes a genePred file using the data (list of lists) given in gene_pred_data

    :param gene_pred_data:          list of lists containing all entries of the genePred file
    :param gene_pred_file_handle:   file handle to the genePred file
    :return:
    """

    print("writing genePred file...")

    for line in gene_pred_data:
        gene_pred_file_handle.write("\t".join(line))

    # close genePred file
    gene_pred_file_handle.close()

    return


def extract_genome_file(genome_file_path, working_directory):
    """
                extracts a IGV .genome file in the working directory

    :param genome_file_path:        file path to the .genome file
    :param working_directory:       path to the working directory
    :return:
    """

    print("extracting genome file...")

    # generate temp folder to extract the files
    temp_folder = os.path.join(working_directory, ".".join(os.path.split(genome_file_path)[1].split('.')[:-1]))
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)

    # extract .genome file
    genome_file = zipfile.ZipFile(genome_file_path, 'r')
    genome_file.extractall(temp_folder)
    genome_file.close()

    return temp_folder


def modify_genome_property_file(property_file_path):
    """
                modifies the property.txt in the extracted genome file by changing the id and name and returning the
                gene file name

    :param property_file_path:  path to the property.txt
    :return:                    file name of the gene file
    """

    print("modifying genome property file (property.txt)...")

    # read property.txt
    with open(property_file_path, 'r') as property_file:
        property_file_content = property_file.readlines()

    # parse content:
    properties = OrderedDict(
        [(line.split('=')[0], "=".join(line.split('=')[1:]).strip()) for line in property_file_content if
         not line.split('=')[0].strip() is ''])

    # modify id and name
    properties["id"] += "_ensembl"
    properties["name"] += " ensembl"

    # write property.txt back to disk
    with open(property_file_path, 'w') as property_file:
        for entry in properties:
            property_file.write(entry + "=" + properties[entry] + "\n")

    return properties["geneFile"]


def replace_gene_file(original_gene_file, generated_gene_file):
    """
                replaces the gene file from the .genome file with the previously generated one

    :param original_gene_file:      file path to the original genePred file (extracted from the .genome file)
    :param generated_gene_file:     file path to the previously generated and modified genePred file
    :return:
    """

    print("replacing gene file...")

    copyfile(generated_gene_file, original_gene_file)


def generate_genome_file(src_folder, output_genome_file_name):
    """
                generates a IGV .genome file from all files in src_folder

    :param src_folder:          path to directory containing all required files
    :param output_genome_file_name:  file path for the output genome
    :return:
    """

    print("generating genome file...")

    # get all files from src_folder
    files = [f for f in os.listdir(src_folder) if os.path.isfile(os.path.join(src_folder, f))]

    # generate .genome (zip) file
    genome_file = zipfile.ZipFile(output_genome_file_name, "w", zipfile.ZIP_DEFLATED)
    for f in files:
        genome_file.write(os.path.join(src_folder, f), f)
    genome_file.close()

    return


def generate_temp_file(named=False, mode='w+t', delete=True):
    """
                generates a (named) temporary file

    :param named:   if true a NamedTemporaryFile instead of a TemporaryFile is created
                        (can be necessary for file system operations)
    :param mode:    defines the mode of the created file (default: 'w+t')
    :return:
    """

    if named:
        temp_file = tempfile.NamedTemporaryFile(mode=mode, delete=delete)
    else:
        temp_file = tempfile.TemporaryFile(mode=mode, delete=delete)

    return temp_file


def main():
    args = parse_args()

    # exit script if input file is missing
    if not validate_args(args):
        return

    ### convert gff file to genePred
    temp_files = {}
    # read input gff
    header, content = read_gff_file(args.gff_file)

    # sort input gff
    sorted_content = sort_gff(content)

    # write sorted gff to file
    temp_files["sorted gff file"] = generate_temp_file(True, 'w+t', False)
    write_gff(temp_files["sorted gff file"], header, sorted_content)

    # setup converter
    temp_files["gff3ToGenePred binary"] = generate_temp_file(True, 'w+b', False)
    setup_gff_converter(temp_files["gff3ToGenePred binary"])

    # run converter
    temp_files["genePred file"] = generate_temp_file(True, 'w+t', False)

    temp_files["genePred file"] = run_gff_converter(temp_files["gff3ToGenePred binary"], temp_files["sorted gff file"],
                                                    temp_files["genePred file"])

    ### modify genePred file
    # parse HGNC file:
    gene_to_hgnc, hgnc_to_gene = read_hgnc_file(args.hgnc_file)

    # generate ENSG-HGNC mapping
    ensg_to_hgnc, hgnc_to_ensg, ensg_to_non_hgnc_gene, non_hgnc_gene_to_ensg = \
        generate_ensg_hgnc_mapping(args.gff_file, hgnc_to_gene.keys())

    # read genePred
    gene_pred_data = read_gene_pred_file(temp_files["genePred file"])

    # modify genePred to fit IGV requirements
    modified_gene_pred_data = modify_gene_pred_data(gene_pred_data, ensg_to_hgnc, hgnc_to_gene, ensg_to_non_hgnc_gene)

    # write modified genePred data back to disk
    temp_files["genePred file modified"] = generate_temp_file(True, 'w+t', False)
    write_gene_pred_file(modified_gene_pred_data, temp_files["genePred file modified"])

    ### replace gene file in genome file
    # generate temporary directory
    temp_files["zip extraction folder"] = tempfile.mkdtemp()
    # extract given .genome file
    temp_folder = extract_genome_file(args.genome_file, temp_files["zip extraction folder"])
    # modify property.txt
    gene_file_name = modify_genome_property_file(os.path.join(temp_folder, "property.txt"))
    # replace gene file
    replace_gene_file(os.path.join(temp_folder, gene_file_name), temp_files["genePred file modified"].name)
    # repack the files into a .genome file
    generate_genome_file(temp_folder, args.output)

    ### cleanup temp files
    # remove generated genePred file
    os.remove(temp_files["genePred file modified"].name)
    # remove zip extraction folder
    shutil.rmtree(temp_files["zip extraction folder"])

    return


if __name__ == '__main__':
    main()
