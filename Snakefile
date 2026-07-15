configfile: "config/config.yaml"
OUT = config["project"]["output_root"]

rule all:
    input: f"{OUT}/localsamap/mapping_statistics.csv", f"{OUT}/localsamap/run_metadata.json"

rule validate_inputs:
    input: "config/config.yaml"
    output: f"{OUT}/validation/input_report.json"
    shell: "python -m workflow.scripts.validate_inputs --config {input} --output {output}"

rule original_samap:
    input: validation=f"{OUT}/validation/input_report.json", config="config/config.yaml"
    output: model=f"{OUT}/original_samap/samap.pkl", scores=f"{OUT}/original_samap/mapping_scores.csv", metadata=f"{OUT}/original_samap/run_metadata.json"
    shell: "python -m workflow.scripts.run_original_samap --config {input.config}"

rule localsamap:
    input: model=f"{OUT}/original_samap/samap.pkl", baseline=f"{OUT}/original_samap/mapping_scores.csv", config="config/config.yaml"
    output: scores=f"{OUT}/localsamap/mapping_statistics.csv", metadata=f"{OUT}/localsamap/run_metadata.json"
    shell: "python -m workflow.scripts.run_localsamap --config {input.config} --model {input.model}"
