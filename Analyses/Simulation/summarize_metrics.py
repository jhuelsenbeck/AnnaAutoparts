import glob

metric_files = glob.glob("metrics_*")

file_summaries = {}

for file_name in metric_files:
    with open(file_name, "r") as f:
        file_name_tokens = file_name.split("_")
        file_name_label = file_name_tokens[1] + "_" + file_name_tokens[2]
        summary = {}
        counts = {}
        lines = f.readlines()
        for l in lines:
            tokens = l.split(" ")
            if "," in l:
                label = tokens[0].split("(")[0]
                if label in counts:
                    summary[label] += float(tokens[3])
                    counts[label] += 1
                else:
                    summary[label] = float(tokens[3])
                    counts[label] = 1
            elif l[0:3] == "RGF" or "(Alpha)" in l:
                summary[tokens[0]] = float(tokens[3])
        for key in counts:
            summary[key] /= counts[key]

        if file_name_label in file_summaries:
            for key in file_summaries[file_name_label]:
                try:
                    file_summaries[file_name_label][key] = (file_summaries[file_name_label][key] + summary[key]) / 2.0
                except:
                    pass
        else:
            file_summaries[file_name_label] = summary

for file in file_summaries:
    print(f"FILE: {file}")
    for key in file_summaries[file]:
        print(f"\t{key}\t\t{file_summaries[file][key]}")
