#!/usr/bin/env cwltool
cwlVersion: v1.2
class: Workflow

requirements:
  - class: StepInputExpressionRequirement
  - class: ScatterFeatureRequirement

inputs:
  ref_genome: File
  query_genes: File[]

outputs:
  blastdb:
    type: Directory
    outputSource: mkdir/blastdb
  blast_results:
    type: File[]
    outputSource: blastn/blast_results

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
            - entry: $(inputs.ref_genome)
              writable: False
      #makeblastdb -in RefGenome -dbtype nucl
      baseCommand: makeblastdb
      inputs:
        ref_genome:
          type: File
          inputBinding:
            prefix: -in
        dbtype:
          type: string?
          default: nucl
          inputBinding:
            prefix: -dbtype
        out_prefix:
          type: string
          inputBinding:
            prefix: -out

      outputs:
        blastfiles:
          type: File[]
          outputBinding:
            glob: "*"
    in:
      ref_genome: ref_genome
      out_prefix:
        source: ref_genome
        valueFrom: $(self.nameroot)
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
            mkdir -p $(inputs.blastdir) && cp
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
            glob: "blast"

    in:
      blastfiles: makeblastdb/blastfiles
      blastdir:
        source: ref_genome
        valueFrom: blast/blastdb/$(self.nameroot)
    out:
      [blastdb]

  blastn:
    run:
      class: CommandLineTool
      # hints:
      #   DockerRequirement:
      #     dockerPull: ncbi/amr:18.06
      requirements:
        - class: InitialWorkDirRequirement
          listing:
            - entry: $(inputs.blastdb)
              writable: False
      # blastn -db -query -out -outfmt
      baseCommand: blastn
      inputs:
        blastdb:
          type: Directory
          inputBinding:
            prefix: -db
            valueFrom: blast/blastdb/$(inputs.ref_genome.nameroot)/$(inputs.ref_genome.nameroot)
        query:
          type: File
          inputBinding:
            prefix: -query
        ref_genome:
          type: File
        outfmt:
          type: int?
          default: 6
          inputBinding:
            prefix: -outfmt
      outputs:
        blast_results:
          type: stdout
      stdout: $(inputs.query.nameroot)_blast.txt
    scatter: query
    in:
      query: query_genes
      blastdb: mkdir/blastdb
      ref_genome: ref_genome
    out:
      [blast_results]
