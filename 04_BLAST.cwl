#!/usr/bin/env cwltool
cwlVersion: v1.2
class: Workflow

requirements:
  - class: StepInputExpressionRequirement

inputs:
  fasta_check_dummy: File?
  gff_check_dummy: File?
  fasta: File
  # query_genes: File[]

outputs:
  blastdb:
    type: Directory
    outputSource: mkdir/blastdb
  # blast_result:
  #   type: File
  #   outputSource: blastn/blast_result

steps:
  makeblastdb:
    run:
      class: CommandLineTool
      # hints:
      #   DockerRequirement:
      #     dockerPull: ncbi/amr:18.06
      requirements:
        - class: InitialWorkDirRequirement
          listing:
            - entry: $(inputs.fasta)
              writable: False
      #makeblastdb -in RefGenome -dbtype nucl
      baseCommand: makeblastdb
      inputs:
        fasta:
          type: File
          inputBinding:
            prefix: -in
        dbtype:
          type: string?
          default: nucl
          inputBinding:
            prefix: -dbtype

      outputs:
        blastfiles:
          type: File[]
          outputBinding:
            glob: "*"
    in:
      fasta: fasta
    out:
      [blastfiles]

  mkdir:
    run:
      class: CommandLineTool
      requirements:
        - class: ShellCommandRequirement
      arguments:
        - shellQuote: false
          valueFrom: >-
            mkdir $(inputs.blastdir) && cp
      inputs:
        blastfiles:
          type: File[]
          inputBinding:
            position: 1
        blastdir:
          type: string
          inputBinding:
            position: 2

      outputs:
        blastdb:
          type: Directory
          outputBinding:
            glob: "$(inputs.blastdir)"

    in:
      blastfiles: makeblastdb/blastfiles
      blastdir:
        source: fasta
        valueFrom: $(self.nameroot)
    out:
      [blastdb]

  # blastn:
  #   run:
  #     class: CommandLineTool
  #     # hints:
  #     #   DockerRequirement:
  #     #     dockerPull: ncbi/amr:18.06
  #   requirements:
  #     - class: InitialWorkDirRequirement
  #       listing:
  #         - entry: $(inputs.fasta)
  #           writable: False
  #   # blastn -db -query -out -outfmt
  #   baseCommand: blastn
  #   inputs:
  #     blastdb:
  #       type: Directory
  #       inputBinding:
  #         prefix: -db
  #     query:
  #       type: File
  #       inputBinding:
  #         prefix: -query
  #     outfmt:
  #       type: int
  #       inputBinding:
  #         prefix: -outfmt

  #   outputs:
  #     blast_result:


