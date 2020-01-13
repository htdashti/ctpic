import os
import main


def clean_folder(entry_path):
    files = ["output.json", "msg", "outputs.zip"]
    for _ in files:
        os.system("rm -rf %s" % os.path.join(entry_path, _))


root_db = "../data/Components-pub_alatis_processed/"

fout_missing_sdf = open("../missing_sdf.txt", "w")
fout_failed_sdf = open("../failed_sdf.txt", "w")
fout_best_penalties = open("../best_penalties.txt", "w")
for an_entry in sorted(os.listdir(root_db)):
    check_id = ""
    if check_id:
        an_entry = check_id
    print(an_entry)
    entry_folder = os.path.join(os.path.abspath(root_db), an_entry)
    clean_folder(entry_folder)
    cmp_name = "alatis_output_compound.sdf"
    sdf_path = os.path.join(entry_folder, cmp_name)
    if not os.path.exists(sdf_path):
        fout_missing_sdf.write("%s\n" % an_entry)
        continue
    best_penalty = main.run_on_a_structure(output_root=entry_folder, sdf_file_path=sdf_path)
    fout_best_penalties.write("%s\t%.02f\n" % (an_entry, best_penalty))
    if check_id:
        break
fout_best_penalties.close()
fout_missing_sdf.close()
fout_failed_sdf.close()
