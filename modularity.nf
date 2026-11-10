nextflow.enable.dsl = 2

params.base_dir   = "mOTUs_gcc"
params.type       = "mb"
params.temps      = "global, low, mid, high"
params.n_null     = 100
params.n_mod_runs = 5
params.outdir     = "try_mod"

process NET_METRICS {
    tag "${temp}"
    publishDir "${params.outdir}/per_temp", mode: 'copy'

    input:
    tuple val(temp), path(edges), path(nodes)

    output:
    path "${params.type}_${temp}_metrics.csv"

    script:
    """
    net_metrics.R --temp ${temp} --edges ${edges} --nodes ${nodes} \\
        --n_null ${params.n_null} --n_mod_runs ${params.n_mod_runs} \\
        --out ${params.type}_${temp}_metrics.csv
    """
}

process BETA_PAIR {
    tag "${t1}_vs_${t2}"
    publishDir "${params.outdir}/beta_pairs", mode: 'copy'

    input:
    tuple val(t1), path(e1, stageAs: 'e1.csv'), path(n1, stageAs: 'n1.csv'),
          val(t2), path(e2, stageAs: 'e2.csv'), path(n2, stageAs: 'n2.csv')

    output:
    path "beta_${t1}_${t2}.csv"

    script:
    """
    beta_pair.R --t1 ${t1} --e1 e1.csv --n1 n1.csv \\
                --t2 ${t2} --e2 e2.csv --n2 n2.csv \\
                --out beta_${t1}_${t2}.csv
    """
}

workflow {
    nets = Channel.fromList(params.temps.tokenize(','))
        .map { t -> tuple(t,
                          file("${params.base_dir}/${t}/edges_${params.type}.csv", checkIfExists: true),
                          file("${params.base_dir}/${t}/nodes_${params.type}.csv", checkIfExists: true)) }

    // una tarea por temperatura (en paralelo)
    NET_METRICS(nets)
        .collectFile(name: "${params.type}_nestedness_modularity.csv",
                     keepHeader: true, skip: 1, storeDir: params.outdir)

    // una tarea por pareja de temperaturas (it[0] < it[3] evita duplicados y auto-parejas)
    pairs = nets.combine(nets).filter { it[0] < it[3] }
    BETA_PAIR(pairs)
        .collectFile(name: "${params.type}_beta_turnover_rewiring.csv",
                     keepHeader: true, skip: 1, storeDir: params.outdir)
}