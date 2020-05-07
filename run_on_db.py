import os
import main
import sys
import json


def clean_folder(entry_path):
    files = ["output.json", "msg", "outputs.zip"]
    for _ in files:
        os.system("rm -rf %s" % os.path.join(entry_path, _))


root_db = os.path.abspath("../data/Components-pub_alatis_processed/")
# root_db = os.path.abspath("/home/nmrbox/hdashti/Desktop/Carbohydrates_2019/codes_for_paper/extract_examples/examples_saccharides")
# root_db = os.path.abspath("/home/nmrbox/hdashti/Desktop/Carbohydrates_2019/codes_for_paper/extract_examples/non_saccharides/")

fout_missing_sdf = open("../missing_sdf.txt", "w")
fout_failed_sdf = open("../failed_sdf.txt", "w")
fout_best_penalties = open("../best_penalties.txt", "w")
for an_entry in sorted(os.listdir(root_db)):
    check_id = ""
    if check_id:
        an_entry = check_id
    print(an_entry)
    entry_folder = os.path.join(root_db, an_entry)
    clean_folder(entry_folder)
    cmp_name = "alatis_output_compound.sdf"
    # cmp_name = "compound.sdf"
    sdf_path = os.path.join(entry_folder, cmp_name)
    if not os.path.exists(sdf_path):
        fout_missing_sdf.write("%s\n" % an_entry)
        continue
    best_penalty = main.run_on_a_structure(output_root=entry_folder, sdf_file_path=sdf_path)
    fout_best_penalties.write("%s\t%.02f\n" % (an_entry, best_penalty))
    """ write msg """
    json_path = os.path.join(entry_folder, "msg.json")
    text = open(os.path.join(entry_folder, "alatis_output_compact.sdf"), "r").read()
    content = text.split("> <ALATIS_Standard_InChI>")
    content = content[1].split("\n")
    inchi = content[0]
    msg = {'sdf_file_path': "alatis_output_compound.sdf",
           "alatis_url": "",
           "alatis_inchi": inchi,
           "input_structure": "alatis_output_compound.sdf",
           "alatis_status": "",
           "error": ""}
    json.dump(msg, open(json_path, "w"), indent=4)

    os.chdir("/home/nmrbox/hdashti/Desktop/Carbohydrates_2019/carbohydrate_identification_git")
    os.system("rm -f obabel_similar_cmps.sdf")
    if check_id:
        break
fout_best_penalties.close()
fout_missing_sdf.close()
fout_failed_sdf.close()

