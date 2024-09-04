source(snakemake@params[['packages']])

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


noise <- snakemake@params[['noise']]

metadata <- fread(snakemake@input[['metadata']])


plot.list <- lapply(noise, function(feature) {
    
    metadata %>%

        ggplot(aes(x=project, y=get(feature))) + 
            geom_violin(fill='steelblue', color='black') + theme_bw() + 
            theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + 
            ylab(feature) + scale_y_log10(labels=scales::label_log())
        
})

plot1 <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=3)

ggsave(plot=plot1, width=15, height=9, filename=snakemake@output[['plot_1']])


plot2 <- metadata %>% 

            ggplot(aes(x=nCount_ATAC, y=nFeature_ATAC)) + geom_point(alpha=0.3) +
            theme_bw() + xlab('Number of UMIs') + ylab('Number of Peaks')


ggsave(plot=plot2, width=15, height=9, filename=snakemake@output[['plot_2']])


g <- metadata %>% 
        ggplot(aes(x=nCount_ATAC)) + geom_histogram(color='black', fill='steelblue', bins=300) + 
        theme_bw() + scale_x_log10(breaks=scales::log_breaks(n=10), label=scales::label_log())


ggsave(plot=g, width=15, height=9, filename=snakemake@output[['plot_3']])


g <- metadata %>% 
        ggplot(aes(x=nFeature_ATAC)) + geom_histogram(color='black', fill='steelblue', bins=200) + 
        theme_bw() + scale_x_log10(breaks=scales::log_breaks(n=10), label=scales::label_log())


ggsave(plot=g, width=15, height=9, filename=snakemake@output[['plot_4']])


g <- metadata %>% 
        ggplot(aes(x=reads_count)) + geom_histogram(color='black', fill='steelblue', bins=300) + 
        theme_bw() + scale_x_log10(breaks=scales::log_breaks(n=10), label=scales::label_log())


ggsave(plot=g, width=15, height=9, filename=snakemake@output[['plot_5']])


g <- metadata %>% 
        ggplot(aes(x=frip)) + geom_histogram(color='black', fill='steelblue', bins=200) + 
        theme_bw() + scale_x_log10(breaks=scales::log_breaks(n=10), label=scales::label_log())


ggsave(plot=g, width=15, height=9, filename=snakemake@output[['plot_6']])







