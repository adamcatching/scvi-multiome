source(snakemake@params[['packages']])

noise <- snakemake@params[['noise']]

metadata <- fread(snakemake@input[['metadata']])[
    predicted_gmm_doublets %chin% 'singlet'
]


plot.list <- lapply(noise, function(feature) {
    
    metadata %>%

        ggplot(aes(x=snakemake@params[['project_name']], y=get(feature))) + 
            geom_violin(fill='steelblue', color='black') + theme_bw() + 
            theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + 
            ylab(feature) + scale_y_log10(labels=scales::label_log())
            
})

p <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=3)

ggsave(plot=p, width=15, height=9, filename=snakemake@output[['plot_1']])


p <- metadata %>% 

        ggplot(aes(x=total_counts, y=n_genes_by_counts)) + geom_point(alpha=0.5) +
        theme_bw() + xlab('Number of UMIs') + ylab('Number of Genes')


ggsave(plot=p, width=15, height=9, filename=snakemake@output[['plot_2']])



g <- metadata %>% 
        ggplot(aes(x=total_counts)) + geom_histogram(color='black', fill='steelblue', bins=300) + 
        theme_bw() + scale_x_log10(breaks=scales::log_breaks(n=10), label=scales::label_log())


ggsave(plot=g, width=15, height=9, filename=snakemake@output[['plot_3']])



g <- metadata %>% 
        ggplot(aes(x=n_genes_by_counts)) + geom_histogram(color='black', fill='steelblue', bins=200) + 
        theme_bw() + scale_x_log10(breaks=scales::log_breaks(n=10), label=scales::label_log())


ggsave(plot=g, width=15, height=9, filename=snakemake@output[['plot_4']])




