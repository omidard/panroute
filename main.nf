#!/usr/bin/env nextflow
/*
 * PanRoute — retrosynthetic route search across all KEGG genomes.
 * Self-contained, parametrised: runs ANY start->end bioconversion.
 *
 *   nextflow run main.nf --queries conf/queries.csv -profile standard
 *   nextflow run main.nf --start C00024 --end C00207 --feedstock C00033
 *
 * The heavy cost (KEGG fetches) is cached on disk under params.cache and reused across
 * runs and queries; use -resume to skip completed queries.
 */
nextflow.enable.dsl = 2

params.queries      = null                 // CSV: name,start,end,feedstock  (feedstock optional)
params.start        = null                 // or a single query via --start/--end
params.end          = null
params.feedstock    = ''
params.name         = 'query'
params.cache        = "$projectDir/cache"
params.assets       = "$projectDir/assets"
params.data         = "/data/bioconversion/data"
params.outdir       = "$projectDir/results"
params.expand_depth = 3
params.max_len      = 12
params.max_routes   = 200
params.min_shared_c = 1
params.no_thermo    = false
params.reactions    = null                 // optional explicit reaction-id file (skips expansion)
params.consensus    = "/data/bioconversion/thermo/directionality_consensus.json"

process PANROUTE {
    tag   "${name}:${start}->${end}"
    publishDir "${params.outdir}/${name}", mode: 'copy'
    maxForks 1                              // serialise to respect KEGG rate limits

    input:
    tuple val(name), val(start), val(end), val(feedstock)

    output:
    tuple val(name), path("report_${start}_${end}.json"), path("routes_${start}_${end}.json"),
          path("species_${start}_${end}.csv"), path("webapp_${start}_${end}.json"),
          path("per_genome_${start}_${end}.csv")

    script:
    def feedArg   = feedstock?.trim() ? "--feedstock ${feedstock}" : ""
    def thermoArg = params.no_thermo ? "--no-thermo" : ""
    def rxnArg    = params.reactions ? "--reactions ${params.reactions}" : ""
    """
    export PYTHONPATH=${projectDir}
    python -m panroute.cli \
        --start ${start} --end ${end} ${feedArg} ${thermoArg} ${rxnArg} \
        --cache ${params.cache} --assets ${params.assets} --data ${params.data} \
        --consensus ${params.consensus} \
        --expand-depth ${params.expand_depth} --max-len ${params.max_len} \
        --max-routes ${params.max_routes} --min-shared-c ${params.min_shared_c} \
        --out .
    """
}

process AGGREGATE {
    tag "aggregate"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path reports

    output:
    path "panroute_summary.json"
    path "panroute_index.json"

    script:
    """
    export PYTHONPATH=${projectDir}
    python ${projectDir}/bin/aggregate.py ${reports} > panroute_summary.json
    python ${projectDir}/bin/build_index.py . > panroute_index.json
    """
}

workflow {
    if (params.queries) {
        queries = Channel.fromPath(params.queries)
            .splitCsv(header: true)
            .map { row -> tuple(row.name, row.start, row.end, row.feedstock ?: '') }
    } else if (params.start && params.end) {
        queries = Channel.of(tuple(params.name, params.start, params.end, params.feedstock))
    } else {
        error "Provide --queries <csv> OR --start <cid> --end <cid>"
    }

    out = PANROUTE(queries)
    AGGREGATE(out.map { it[1] }.collect())
}
