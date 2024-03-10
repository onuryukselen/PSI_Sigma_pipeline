#!/usr/bin/env python3

import sys, os, re, glob, argparse, subprocess
from collections import defaultdict
from textwrap import dedent
from multiprocessing import Pool

def main(args):

	groups = read_groups(args.all_groups, args.compare_file)

	p = Pool(args.threads)

	p.map(filter_gct, ((args.script_path, group, args.min_control, args.min_treatment) for group in groups))
	p.close()
	p.join()

	prepare_barchart_summary(groups)

	p = Pool(args.threads)
	p.map(build_report, ((group, groups[group][0], groups[group][1]) for group in groups))
	p.close()
	p.join()

def read_groups(groups_file, compare_file):
	
	groups = {}

	if compare_file[:7] == 'NO_FILE':
		controls = {}
		treatments = {}
		with open(groups_file) as infile:
			for line in infile:
				sample, code, name = line.split('\t')
				if ''.join(filter(str.isdigit, code)) != '0':
					treatments[name.rstrip()] = ''.join(filter(str.isalpha, code))
				else:
					controls[''.join(filter(str.isalpha, code))] = name.rstrip()

		for treatment in treatments:
			groups[treatment] = (controls[treatments[treatment]], treatment)

	else:
		with open(compare_file) as infile:
			infile.readline()
			for line in infile:
				cur = line.rstrip().split()
				control = cur[0]
				treatment = cur[1]
				name = cur[2]
				groups[name] = (control, treatment)

	return groups

def filter_gct(args):

	script_path, group, min_control, min_treatment = args

	if len(glob.glob('%s*sorted.annotated.txt' % group)) == 1:
		sorted_file = glob.glob('%s*sorted.annotated.txt' % group)[0]
	else:
		sorted_file = glob.glob('%s*sorted.txt' % group)[0]

	gct_file = glob.glob('%s*denominator.gct' % group)[0]
	cmd = "%s %s %s %d %d" % (script_path, sorted_file, gct_file, min_control, min_treatment)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print('Filtering: %s' % group)
	print(out.decode())

def build_report(args):
	
	name, control, treatment = args

	output = []

	volcano_file = glob.glob('%s*.volcano.txt' % name)[0]

	output.append(print_header(name))
	output.append(print_libraries())
	output.append(print_functions())
	output.append(read_volcano_data(volcano_file))
	output.append(volcano_results_table(name))
	output.append(volcano_plot())
	output.append(jitter_plot(control, treatment))
	output.append(read_barchart_data())
	output.append(barchart_results_table())
	output.append(barchart())
	write_output(name, output)
	run_markdown(name)

def print_header(name):

	return(dedent(
	'''
	---
	title: {}
	date: "`r Sys.Date()`"
	output:
	  html_document:
	    code_folding: hide
	---
	'''.format(name)).strip())

def print_libraries():
	
	return(dedent(
	'''
	```{r, load_libraries, message=FALSE, include=FALSE}
	# Load Libraries
	library(ggplot2)
	library(ggrepel)
	library(dplyr)
	library(DT)
	library(scales)
	library(tidyr)
	```'''))

def print_functions():

	return(dedent(
	'''
	```{r, functions, include=FALSE}
	reverselog = function() {
	  trans_new("reverselog", function(x) -log10(x), function(x) 10^(-x), log_breaks(base = 10), domain = c(1e-1000, Inf))
	}
	
	volcano_plot = function(df, padj_cutoff=.01, dpsi_cutoff=20,
	                        positive_color='firebrick', noChange_color='grey', negative_color='steelblue',
	                        positive_alpha=1, noChange_alpha=.3,negative_alpha=1,
	                        positive_size=2, noChange_size=1, negative_size=2,
	                        dpsi_markers=TRUE, dpsi_marker_color='grey', center_marker=TRUE, center_marker_color='black',
	                        padj_marker=TRUE, padj_marker_color='grey',
	                        display_positive_count = TRUE, display_negative_count = TRUE, display_noChange_count=FALSE) {
	
	  colors  = c('Negative'=negative_color, 'No Change'=noChange_color, 'Positive'=positive_color)
	  sizes   = c('Negative'=negative_size,  'No Change'=noChange_size,  'Positive'=positive_size)
	  alphas  = c('Negative'=negative_alpha, 'No Change'=noChange_alpha, 'Positive'=positive_alpha)
	
	  positive_count = nrow(df %>% filter(Group == 'Positive'))
	  noChange_count = nrow(df %>% filter(Group == 'No Change'))
	  negative_count = nrow(df %>% filter(Group == 'Negative'))
	
	  return(
	    ggplot(data, aes(x=dPSI, y=pvalue, color=Group, alpha=Group, size=Group, label=Label)) +
	      theme_classic() +
	      theme(legend.position = 'none') +
	      scale_x_continuous(limits=c(-100, 100), name='dPSI') +
	      scale_y_continuous(trans=reverselog(), name='Significance', labels=trans_format('log10',math_format(10^.x))) +
	      scale_color_manual(values=colors) +
	      scale_alpha_manual(values=alphas) +
	      scale_size_manual(values=sizes) +
	      geom_hline(yintercept=.01, linetype=2, color='grey') +
	      {if (padj_marker) geom_hline(yintercept = padj_cutoff, linetype=2, color=padj_marker_color)} +
	      {if (dpsi_markers) geom_vline(xintercept = dpsi_cutoff, linetype=2, color=dpsi_marker_color)} +
	      {if (center_marker) geom_vline(xintercept = 0, linetype=2, color=center_marker_color)} +
	      {if (dpsi_markers) geom_vline(xintercept = -dpsi_cutoff, linetype=2, color=dpsi_marker_color)} +
	      geom_point() +
	      {if (display_positive_count) annotate('text', x=Inf, y=0, hjust=1, vjust=1, color=positive_color, label=paste0("Positive: ", positive_count))} +
	      {if (display_negative_count) annotate('text', x=-Inf, y=0, hjust=-.01, vjust=1, color=negative_color, label=paste0("Negative: ", negative_count))} +
	      {if (display_noChange_count) annotate('text', x=0, y=0, vjust=1, color=noChange_color, label=paste0("No Change: ", noChange_count))} +
	      geom_label_repel(label.size=NA, fill=NA, na.rm=TRUE, max.overlaps = 50, max.time = 5)
	  )
	}

	jitter_plot = function(data, control, treatment, padj_cutoff=.01,
	                       CE_color='#479FF8', A5SS_color='#EFBD40', A3SS_color='#EA4025', IR_color='#B45084', non_sig_color='grey',
	                       sig_size=.5, non_sig_size=.3,
	                       sig_alpha=1, non_sig_alpha=.3) {
	
	  jitter_data = data %>% 
	    mutate(Event.Type = case_when(Event.Type == 'Exon Inclusion' ~ "CE",
	                                                 Event.Type == 'Exon Skipping' ~ "CE",
	                                                 Event.Type == 'Increased IR' ~ 'IR',
	                                                 Event.Type == 'Decreased IR' ~'IR',
	                                                 Event.Type == "Alt. 3'-splice-site" ~ 'A3SS',
	                                                 Event.Type == "Alt. 5'-splice-site" ~ 'A5SS'
	    )) %>%
	    mutate(Event.Type = factor(Event.Type, levels=c('IR', "A3SS", "A5SS", "CE"))) %>%
	    mutate(Significant = pvalue < padj_cutoff) %>%
	    mutate(Group = paste(Event.Type, Significant, sep='_'))
	
	  colors = c("CE_TRUE"=CE_color, "CE_FALSE"=non_sig_color, "A3SS_TRUE"=A3SS_color, "A3SS_FALSE"=non_sig_color, "A5SS_TRUE"=A5SS_color, "A5SS_FALSE"=non_sig_color, "IR_TRUE"='#B45084', "IR_FALSE"=non_sig_color)
	
	  return(
	    ggplot(jitter_data, aes(x=Event.Type, y=dPSI, alpha=Significant, color=Group, size=Significant)) +
	      theme_classic(base_size = 16) +
	      theme(axis.title.y = element_blank(),
	            axis.ticks.y = element_blank(),
	            axis.line.y = element_blank(),
	            legend.position = 'none') +
	      scale_y_continuous(limits=c(-100,100), name=paste0('∆PSI\n[',treatment,' - ', control, ']')) +
	      scale_color_manual(values=colors) +
	      scale_alpha_manual(values=c(non_sig_alpha, sig_alpha)) +
	      scale_size_manual(values=c(non_sig_size, sig_size)) +
	      geom_jitter(data = jitter_data %>% filter(Significant==FALSE)) +
	      geom_jitter(data = jitter_data %>% filter(Significant==TRUE)) +
	      coord_flip() +
	      geom_boxplot(data=jitter_data, mapping=aes(x=Event.Type, y=dPSI), color='black', fill=NA, inherit.aes = FALSE, outlier.colour = NA)
	  )
	}

	barchart_plot = function(data, control, treatment,
	                        exon_inclusion_color='#479FF8', exon_skipping_color='#81D653',
	                        alt5_color='#EFBD40', alt3_color='#EA4025',
	                        increase_ir_color='#B45084', decrease_ir_color='#5F5F5F') {
	                        
	  df = data %>%
	       pivot_longer(!`File Name`, names_to = "EventType", values_to = "Count") %>%
	       mutate(EventType = factor(EventType, levels=c("Exon Inclusion", "Exon Skipping", "Alt. 5'-splice-site", "Alt. 3'-splice-site", "Increased IR", "Decreased IR")))
	  
	  colors  = c('Exon Inclusion'=exon_inclusion_color, 'Exon Skipping'=exon_skipping_color, "Alt. 5'-splice-site"=alt5_color, "Alt. 3'-splice-site"=alt3_color, 'Increased IR'=increase_ir_color, 'Decreased IR'=decrease_ir_color)
	  
	  return(
	    ggplot(df, aes(x=EventType, y=Count, fill=EventType)) +
	      facet_wrap(~`File Name`, nrow=1, strip.position = 'bottom') +
	      theme_classic() +
	      theme(axis.text.x=element_blank(),
	            axis.ticks.x=element_blank(),
	            axis.title.x=element_blank(),
	            strip.background=element_blank(),
	            strip.placement = 'outside') +
	      scale_y_continuous(expand=c(0,0), name='Number of Genes (|dPSI| > 20 and p-value <0.01)') +
	      scale_fill_manual(values=colors, name='') +
	      geom_bar(stat='identity')
	  )
	}
	```
	'''
	))

def read_volcano_data(file):

	return(dedent(
	'''
	```{{r, read_volcano_data}}
	data = read.delim("{}", header=TRUE, sep='\\t') %>%
	       mutate(pvalue = 10^(-log10.p.value.)) %>%
	       mutate(Significant = case_when(dPSI > 20 & pvalue < .01 ~ 'Significant',
	                                      dPSI < -20 & pvalue < .01 ~ 'Significant',
	                                      TRUE ~ 'Non-significant'
	             )) %>%
	       mutate(Group = case_when(Significant == 'Significant' & dPSI< -20 ~ 'Negative',
	                                Significant == 'Significant' & dPSI > 20 ~ 'Positive',
	                                TRUE ~ 'No Change'
	             )) %>%
	       mutate(Label = case_when(Group == 'Positive' | Group == 'Negative' ~ Gene.Symbol))
	```'''.format(file)))

def volcano_results_table(name):
	
	return(dedent(
	'''
	```{{r volcano_table, warning=FALSE}}
	rank = data %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(-abs(dPSI))) %>% select(Gene.Symbol, Event.Region, Target.Exon, Event.Type, dPSI, pvalue, Database.ID, Rank)
	
	df = data %>% 
	     select(Gene.Symbol, Event.Region, Target.Exon, Event.Type, Database.ID, dPSI, pvalue, Group, Significant) %>%
	     left_join(rank, by=c('Gene.Symbol', 'Event.Region', 'Target.Exon', 'Event.Type', 'Database.ID', 'dPSI', 'pvalue')) %>% 
	     mutate(abs_dpsi=-abs(dPSI)) %>% 
	     arrange(Rank, abs_dpsi) %>%
	     select(Gene.Symbol, Event.Region, Target.Exon, Event.Type, dPSI, pvalue, Group)
	
	datatable(df,
	  rownames=FALSE,
	  colnames = c("Gene", "Event Region", "Target Exon", "Event Type", "dPSI", "p-value", "Group"),
	  extensions = 'Buttons',
	  options=list(
	    columnDefs=list(list(visible=FALSE, targets=c(6))),
	    dom = 'lftBipr',
	    buttons = list(
	      list(extend = 'csvHtml5', text='Download', filename = paste0('{}', "_PSI-Sigma_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\\t')
	     )
	    ),
	) %>% 
	formatSignif(columns=c("pvalue"), digits=4) %>%
	formatStyle('Group',target = 'row', color = styleEqual(c("No Change", "Positive", "Negative"), c('black', 'firebrick', 'steelblue')))
	```
	'''.format(name)
	))

def volcano_plot():
	
	return(dedent(
	'''
	```{r, volcano_plot, warning=FALSE}
	volcano_plot(data)
	```
	'''
	))

def jitter_plot(control, treatment):
	
	return(dedent(
	'''
	```{{r, jitter_plot, warning=FALSE}}
	jitter_plot(data, "{}", "{}")
	```
	'''.format(control, treatment)
	))

def read_barchart_data():

	return(dedent(
	'''
	```{r, read_barchart_data, warning=FALSE}
	barchart_data = read.delim('barchart_summary.tsv', sep='\\t', header=TRUE, check.names = FALSE)
	```
	'''
	))

def barchart_results_table():

	return(dedent(
	'''
	```{r, barchart_table}
	datatable(barchart_data,
	  rownames=FALSE,
	  extensions = 'Buttons',
	  options=list(
	    dom = 'lftBipr',
	    buttons = list(
	      list(extend = 'csvHtml5', text='Download', filename = 'barchart_data.tsv', extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	     )
	    )
	)
	```
	'''
	))

def barchart():

	return(dedent(
	'''
	```{r, barchart}
	barchart_plot(barchart_data)
	```
	'''
	))

def prepare_barchart_summary(groups):

	barchart_output = []

	for group in groups:
		
		barchart_file = glob.glob('%s*barchart.txt' % group)[0]
		
		with open(barchart_file) as infile:
			header = infile.readline()
			content = infile.readline().rstrip().split('\t')
			barchart_output.append('%s\t%s' % (group, '\t'.join([i for i in content[1:]])))

	with open('barchart_summary.tsv', 'w') as out:
		out.write(header)
		out.write('\n'.join(barchart_output))

def write_output(name, output):

	with open('%s.Rmd' % (name), 'w') as out:
		out.write('\n'.join(output)) 

def run_markdown(name):

	cmd = "Rscript -e 'rmarkdown::render(\"%s.Rmd\", \"html_document\")'" % (name)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print('%s:' % name)
	print(out.decode())
	print('%s:' % name, file=sys.stderr)
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="summarise_PSI-Sigma.py", description='', usage='%(prog)s [options]')
	
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-a', '--all-groups', required=True, help='Name of all_groups file.', metavar='', dest='all_groups')
	input_args.add_argument('-c', '--comparison-file', required=True, help='Name of comparison file.', metavar='', dest='compare_file')
	input_args.add_argument('-s', '--script-path', default = 'PSI-Sigma_filter_v1.2.pl', help='Path to PSI-Sigma_filter script', metavar='', dest='script_path')
	input_args.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for knitting.', metavar='', dest='threads')

	filter_args = parser.add_argument_group('Filtering')
	filter_args.add_argument('-m', '--min-control', default=2, type=int, help='Minimal number of control samples.', metavar='', dest='min_control')
	filter_args.add_argument('-n', '--min-treatment', default=2, type=int, help='Minimal number of treatment samples.', metavar='', dest='min_treatment')

	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)