#!/usr/bin/env Rscript
# PanRoute figures for: 2,3-butanediol from acetate.
#  A) overview funnel (terminal-enzyme -> encodes-route -> +feedstock) + composition + producers
#  B) all native routes pyruvate -> 2,3-butanediol (carbon-skeleton DAG, enzyme-labelled)
# House standard: Wong palette, 10 pt font floor, biological gene names, white, editable SVG.
suppressMessages({library(ggplot2);library(dplyr);library(patchwork);library(ggtext)
  library(igraph);library(ggraph);library(tidygraph);library(scales);library(jsonlite)})

D  <- "/data/bioconversion/panroute/figs/data"
OUT<- "/data/bioconversion/panroute/figs"
INK<-"#1B2631"; GREY<-"#566573"; PAPER<-"#FFFFFF"
GPOS<-"#0072B2"; GNEG<-"#D55E00"; ARCH<-"#009E73"; OTH<-"#9AA6B2"
ACCENT<-"#117A65"; T0C<-"#AEB6BD"; T2C<-"#5AA9C9"; T3C<-"#117A65"
NB <- " "   # plain space (this gridtext build does not decode &nbsp;)

meta   <- fromJSON(file.path(D,"C03044_meta.json"))
funnel <- read.delim(file.path(D,"C03044_funnel.tsv"))
gram   <- read.delim(file.path(D,"C03044_gram.tsv"))
notable<- read.delim(file.path(D,"C03044_notable.tsv"))
nodes  <- read.delim(file.path(D,"C03044_nodes.tsv"))
edges  <- read.delim(file.path(D,"C03044_edges.tsv"))

base <- theme_minimal(base_size=12)+theme(
  text=element_text(color=INK),
  plot.title=element_markdown(size=14,face="bold",color=INK),
  plot.subtitle=element_markdown(size=11,color=GREY),
  axis.text=element_text(size=11,color=INK), axis.title=element_text(size=11),
  panel.grid.minor=element_blank(),
  plot.background=element_rect(fill=PAPER,color=NA),
  panel.grid.major.y=element_blank())

## ============================ FIGURE A — OVERVIEW =========================
funnel$tier <- factor(funnel$tier, levels=c("T3","T2","T0"))       # T0 top
ylabs <- c(T0="terminal enzyme\n(last step only)",
           T2="encodes a full\nnative route",
           T3="+ can take\nup acetate")
funnel$col  <- c(T0=T0C,T2=T2C,T3=T3C)[as.character(funnel$tier)]
pA <- ggplot(funnel, aes(species, tier, fill=I(col)))+
  geom_col(width=.6)+
  geom_text(aes(label=comma(species)), hjust=-0.12, color=INK, fontface="bold", size=4.8)+
  scale_y_discrete(labels=ylabs)+
  scale_x_continuous(expand=expansion(mult=c(0,.16)))+
  labs(title=paste0("**a**",NB,"How many KEGG prokaryotes could make 2,3-butanediol from acetate"),
       subtitle=paste0("honest funnel — gene *presence*, not proof of production",NB,"·",NB,
                       "**",meta$overflow_excluded,"** overflow-only species excluded"),
       x="prokaryote species (strain-collapsed)", y=NULL)+
  base+theme(axis.text.y=element_text(size=11,color=INK,hjust=1,lineheight=.95),
             plot.margin=margin(6,14,6,6))

# gram composition of T2 (donut)
gram$group <- factor(gram$group, levels=c("Gram-positive","Gram-negative","Archaea","Other"))
gcol <- c("Gram-positive"=GPOS,"Gram-negative"=GNEG,"Archaea"=ARCH,"Other"=OTH)
pB <- ggplot(gram, aes(x=2, y=species, fill=group))+
  geom_col(width=1, color="white")+ coord_polar(theta="y")+ xlim(.5,2.5)+
  scale_fill_manual(values=gcol, name=NULL)+
  geom_text(aes(label=species), position=position_stack(vjust=.5), color="white", fontface="bold", size=4.2)+
  labs(title=paste0("**b**",NB,"Composition (n=",comma(meta$T2),")"))+
  theme_void(base_size=12)+
  theme(plot.title=element_markdown(size=13,face="bold"),
        legend.text=element_text(size=11), legend.position="bottom",
        legend.key.size=unit(5,"mm"),
        plot.background=element_rect(fill=PAPER,color=NA),
        plot.margin=margin(2,4,2,4))

# notable species lollipop (validated producers flagged)
notable$species <- factor(notable$species, levels=rev(notable$species))
notable$gcol <- gcol[ifelse(notable$gram=="Gpos","Gram-positive",
                     ifelse(notable$gram=="Gneg","Gram-negative",
                     ifelse(notable$gram=="Arch","Archaea","Other")))]
notable$val <- ifelse(notable$flag=="validated","✓ validated","")
pC <- ggplot(notable, aes(x=1, y=species))+
  geom_point(aes(color=I(gcol)), size=3.6)+
  geom_text(aes(x=1.08, label=species), hjust=0, fontface="italic", size=3.7, color=INK)+
  geom_text(aes(x=3.15, label=val), hjust=0, size=3.6, color=ACCENT, fontface="bold")+
  scale_x_continuous(limits=c(.9,3.9))+
  labs(title=paste0("**c**",NB,"Example route-carrying species"))+
  theme_void(base_size=12)+
  theme(plot.title=element_markdown(size=13,face="bold"),
        plot.background=element_rect(fill=PAPER,color=NA),
        plot.margin=margin(2,6,2,6))

facts <- paste0(
  "<span style='font-size:11pt'>",
  "**Feedstock** acetate → pyruvate",NB,"·",NB,
  "**", meta$n_routes, " routes** (min ", meta$shortest, " steps)",NB,"·",NB,
  "**Validation** P ", sprintf('%.2f',meta$precision), " / R ", sprintf('%.2f',meta$recall),
  "</span>")
pFacts <- ggplot()+annotate("richtext", x=0,y=0,label=facts, hjust=0, vjust=.5,
                            fill="#EEF3F1", label.color=NA, label.padding=unit(c(5,8,5,8),"pt"))+
  xlim(0,1)+ylim(-1,1)+theme_void()+theme(plot.background=element_rect(fill=PAPER,color=NA))

capA <- paste0("**Genome *potential*, not production.** A gene hit means the enzyme is encoded — not expressed, ",
  "regulated, or carrying flux. KEGG genomes are culture-biased; Gram/taxonomy is heuristic. ",
  "The terminal-enzyme count (", comma(meta$T0), ") overcounts because it ignores whether the cell can reach the precursor.")

figA <- (pA / (pB | pC) / pFacts) + plot_layout(heights=c(1.05,1.3,.16)) +
  plot_annotation(
    title="2,3-Butanediol from acetate: native capacity across KEGG prokaryote genomes",
    subtitle="PanRoute retrosynthetic search · carbon-skeleton network · thermodynamically screened · direction-aware feedstock gating",
    caption=capA,
    theme=theme(plot.title=element_markdown(size=15,face="bold",color=INK),
                plot.subtitle=element_markdown(size=11,color=GREY),
                plot.caption=element_textbox_simple(size=10,color=GREY,margin=margin(t=8)),
                plot.background=element_rect(fill=PAPER,color=NA),
                plot.margin=margin(12,16,10,14)))

ggsave(file.path(OUT,"fig_butanediol_overview.pdf"), figA, width=250, height=215, units="mm", device=cairo_pdf)
ggsave(file.path(OUT,"fig_butanediol_overview.png"), figA, width=250, height=215, units="mm", dpi=300, bg="white")
suppressMessages(ggsave(file.path(OUT,"fig_butanediol_overview.svg"), figA, width=250, height=215, units="mm"))
message("OK figA")

## ============================ FIGURE B — PATHWAYS =========================
# canonical = the textbook acetolactate route alsS -> alsD -> budC (not merely any 3-step path)
CANON_ENZ <- c("alsS","alsD","budC")
edges$kind <- ifelse(edges$enzymes %in% CANON_ENZ | grepl("alsS|alsD|budC", edges$enzymes),
                     "canonical: alsS → alsD → budC", "alternate route")
# label only canonical edges (reduce label count, don't shrink); alternates stay grey
edges$elab <- ifelse(edges$kind=="canonical: alsS → alsD → budC", edges$enzymes, "")
g <- graph_from_data_frame(d=edges[,c("from","to","enzymes","elab","reaction","kind")],
                           vertices=nodes[,c("cid","name","role","min_len")], directed=TRUE)
lay <- create_layout(g, layout="sugiyama")
lay2 <- lay; lay2$x <- lay$y; lay2$y <- -lay$x        # flow left -> right

rolecol <- c(start="#1B2631", end=ACCENT, intermediate="#FFFFFF")
pPath <- ggraph(lay2)+
  geom_edge_fan(aes(edge_color=kind, edge_width=kind, label=elab),
                angle_calc="along", label_dodge=unit(3,"mm"),
                label_size=4.0, family="sans", fontface="italic",
                arrow=arrow(length=unit(2.4,"mm"), type="closed"),
                end_cap=circle(7,"mm"), start_cap=circle(7,"mm"),
                label_colour=ACCENT)+
  geom_node_point(aes(fill=role), shape=21, size=8.5, stroke=.7, color=INK)+
  geom_node_label(aes(label=name), size=3.7, repel=TRUE, color=INK, fontface="bold",
                  fill="white", label.size=0, label.padding=unit(.6,"mm"),
                  box.padding=.9, point.padding=.5, seed=7, max.overlaps=Inf)+
  scale_edge_color_manual(values=c("canonical: alsS → alsD → budC"=ACCENT,"alternate route"="#AEB8C0"), name=NULL)+
  scale_edge_width_manual(values=c("canonical: alsS → alsD → budC"=1.2,"alternate route"=.5), name=NULL)+
  scale_fill_manual(values=rolecol, guide="none")+
  labs(title="All native routes from pyruvate to 2,3-butanediol found by PanRoute",
       subtitle=paste0("carbon-skeleton retrosynthetic search over KEGG",NB,"·",NB,meta$n_routes,
                       " routes",NB,"·",NB,"stereoisomers merged",NB,"·",NB,
                       "**canonical route highlighted**"),
       caption="Nodes = metabolites; directed edges = enzymatic steps (biological gene symbols). All routes converge on acetoin, the committed precursor of 2,3-butanediol (final step budC). Alternate routes are other carbon-skeleton paths the search found to reach that convergence.")+
  theme_void(base_size=12)+
  theme(plot.title=element_markdown(size=15,face="bold",color=INK),
        plot.subtitle=element_markdown(size=11,color=GREY),
        plot.caption=element_textbox_simple(size=10,color=GREY,margin=margin(t=6)),
        legend.text=element_text(size=11), legend.position="top",
        plot.background=element_rect(fill=PAPER,color=NA),
        plot.margin=margin(12,18,10,18))

ggsave(file.path(OUT,"fig_butanediol_pathways.pdf"), pPath, width=300, height=200, units="mm", device=cairo_pdf)
ggsave(file.path(OUT,"fig_butanediol_pathways.png"), pPath, width=300, height=200, units="mm", dpi=300, bg="white")
suppressMessages(ggsave(file.path(OUT,"fig_butanediol_pathways.svg"), pPath, width=300, height=200, units="mm"))
message("OK figB")
