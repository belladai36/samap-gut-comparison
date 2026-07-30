configfile: "config/config.yaml"

OUT = config["project"]["output_root"]


rule all:
    input:
        f"{OUT}/original_samap/samap.pkl",
        f"{OUT}/original_samap/mapping_scores.csv",
        f"{OUT}/original_samap/run_metadata.json",


rule validate_inputs:
    input:
        "config/config.yaml"
    output:
        f"{OUT}/validation/input_report.json"
    shell:
        "python -m workflow.scripts.validate_inputs --config {input} --output {output}"


rule original_samap:
    input:
        validation=f"{OUT}/validation/input_report.json",
        config="config/config.yaml"
    output:
        model=f"{OUT}/original_samap/samap.pkl",
        scores=f"{OUT}/original_samap/mapping_scores.csv",
        metadata=f"{OUT}/original_samap/run_metadata.json"
    shell:
        "python -m workflow.scripts.run_original_samap --config {input.config}"


# A local-pattern extension should be maintained and installed as a separate
# package before it is added to this public, official-SAMap workflow.
