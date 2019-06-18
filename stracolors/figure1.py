#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 11:15:19 2018

@author: viola


Code to plot the figures in the Figure 1 panel.
Contains: plot_Heritability (baghera vs ldsc),
          plot_Enrichment (enrichment of significant genes),
          plot_HeritabilityPrevalence (prevalence vs number of significant gene and mi vs number of CHG),
          plot_minSNP (minSNP vs BAGHERA),
          plot_significant (h2 and number of significant genes)
In total there are 6 figures saved (both jpeg and pdf)
"""
import matplotlib
matplotlib.use("TkAgg")

from adjustText import adjust_text
import numpy as np
import pandas as pd
#import gseapy as gp
import argparse
import scipy.stats as stats
import seaborn as sns
from matplotlib import pyplot as plt
import statsmodels.stats.multitest as smm
import glob as g
import matplotlib.gridspec as gridspec
from palettable.cmocean.sequential import *

from palettable.colorbrewer.diverging import *
from palettable.colorbrewer.sequential import *
from palettable.cartocolors.sequential import *
from palettable.cubehelix import Cubehelix
import matplotlib.font_manager as font_manager

from matplotlib import rcParams
import scipy.stats


#from bokeh.io import show, output_file
#from bokeh.models import ColumnDataSource
#from bokeh.plotting import figure
#from bokeh.sampledata.sprint import sprint
#from bokeh.sampledata.commits import data
#from bokeh.transform import jitter



# Color settings
color_blu=['#0571b0']
color_red=['#ca0020']
#palette_binary = ['#ca0020','#0571b0']
palette_binary = RdBu_5.mpl_colors[0::4]
print(palette_binary)
palette_sequential_blu = PuBu_8.mpl_colors #PuBu_9
palette_sequential_red = Amp_20.mpl_colors #Reds_9




def plot_Heritability(tab_heritability_file, figs_output):
    
    """
    Plots the comparison between BAGHERA and LDSC estimates of heritability
    """
    
    with open(tab_heritability_file,"r") as f:
        tab_heritability=pd.read_table(f,sep=",")
    
    tab_heritability=tab_heritability.sort_values(by=["mi_median"], ascending=False)
    
    df=pd.DataFrame(columns=["Malignancy","method","h2_observed"])
    for i,t in tab_heritability.iterrows():
        
        df=df.append({"Malignancy":t["name"],"method":"LDSC","h2_observed":t["h2_observed"]}, ignore_index=True)
        df=df.append({"Malignancy":t["name"],"method":"BAGHERA","h2_observed":t["mi_median"]}, ignore_index=True)
    

    
    #fig,axes=plt.subplots(1, 1, figsize=(10, 10),sharey=True)
    #fig.subplots_adjust(left=0.2,right=0.99) 
    #pal=sns.cubehelix_palette(100)
       
    g = sns.catplot(x="Malignancy", y="h2_observed", hue="method", data=df,
                    height=6, kind="bar", palette=palette_binary, aspect=2)

    g.set_ylabels("Observed Heritability")

    #g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 6)
    g.set_xticklabels(rotation = 90, fontsize = 8) 

    g.despine(bottom=True,left=True)
    

    #g.set_title("Actionable targets")
    g.savefig(figs_output+"ldsc_vs_baghera.pdf", format='pdf')
    
    g.savefig(figs_output+"ldsc_vs_baghera.png", format='png')


def plot_Enrichment(group, figs_output):
    
    # TODO: plot enrichment
    """ Distribution of significant bg for each cancer. Plot saved in the selected folder.
    """
    
    # Initialize the FacetGrid object

    #group=group[group["mi_median"]>0.001].copy()
    group2=group[group["P"]>0.99].copy()
    group2=group2.sort_values(by=['weight'])
    group_weight=group2[group2['weight']>170].sort_values(by=['weight'])
    perc_90=np.percentile(group2['weight'], 90)

    name_temp=[]
    weight_temp=[]
    for i, k in group2.groupby(by=['name']):
        name_temp.append(i)
        weight_temp.append(np.max(k['weight'].values))
    df_temp=pd.DataFrame()
    df_temp['name']=name_temp
    df_temp['weight']=weight_temp
    print('n genes with max enrichment >10 ', len(df_temp[df_temp['weight']>10]))
    print('n genes tot ', len(df_temp))


    print('n genes with enrichment >10 ', len(group2[group2['weight']>10]))
    print('n genes  ', len(group2))

    el_perc_90=group2.iloc[0,:]

    f,axes=plt.subplots(1, figsize=(10,5))
    a=sns.distplot(group2['weight'],color=palette_binary[1] , ax=axes)
    axes.scatter(group_weight['weight'],[0]*len(group_weight), color=palette_binary[0])


    texts=[]
    
    #texts.append(plt.text(el_perc_90['weight'], 0, 'Percentile'+el_perc_90['name'].replace(';','\n')+'\n'+el_perc_90['Malignancy'], fontsize=6 , rotation=-330) )


    for key, row in group_weight.iterrows():
        if row['weight']>=170:
            k=row['name'].split(';')
            texts.append(plt.text(row['weight'], 0, row['name'].replace(';','\n')+'\n'+row['Malignancy'], fontsize=10 , rotation=-330) )
            #axes.annotate(row['name']+'-'+row['Malignancy'], xy=(row['weight'], 0),
            #            xytext=(row['weight'], 0), 
            #            fontsize=6,rotation=-330)   

    adjust_text(texts, ha='center', va='bottom', only_move={'points':'y', 'text':'y'},arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3", lw=0.2))

    axes.scatter(el_perc_90['weight'],[0], color=palette_binary[0])
    axes.annotate(el_perc_90['name'].replace(';','\n')+'\n'+el_perc_90['Malignancy'], xy=(el_perc_90['weight'], 0),
                           xytext=(100, 0.01), 
                            fontsize=10,rotation=-330,arrowprops=dict(arrowstyle="->",
                           connectionstyle="arc3", lw=0.2))   

    axes.set_ylabel("Density", fontsize=10)
    axes.set_xlabel("Fold Change",fontsize=10)
    #axes.set(ylabel="", xlabel='Fold Change')


    f.savefig(figs_output+"CHGs_enrichment.pdf", format='pdf')
    f.savefig(figs_output+"CHGs_enrichment.png", format='png')   

    group2["Enrichment"]=group2["weight"]
    group2["hue"]=1*(group2["mi_median"]>0.01)+1*(group2["mi_median"]>0.001)
    
    malignancies=[]
    mi=[]
    bg_min=[]
    bg_max=[]
    for mal,tab in group2.groupby(["Malignancy"]):
        malignancies.append(mal)
        mi.append(tab["mi_median"].values[0])
        bg_min.append(np.min(tab["bg_median"]))
        bg_max.append(np.max(tab["bg_median"]))
        
    df=pd.DataFrame()
    df["mal"]=malignancies
    df["mi"]=mi
    df["bg_min"]=bg_min
    df["bg_max"]=bg_max
    df=df.sort_values(by=["mi"],ascending=True)
    
    print(df)

    f, ax = plt.subplots(figsize=(7, 6))   

    first_df=df.iloc[-6:,:]
    first_group=group2[group2['Malignancy'].isin(first_df['mal'].values.tolist())]    

    sns.violinplot(x="weight", y="Malignancy", data=first_group, order=first_df['mal'].values.tolist(),inner='points',
                color=palette_sequential_blu[-3],size=2)


    #sns.stripplot(x="weight", y="Malignancy", data=first_group, order=first_df['mal'].values.tolist(),
    #            color='red')
    

    #sns.boxplot(x="bg_median", y="Malignancy", data=group2,order=df["mal"].values.tolist(),
    #                   color=palette_sequential_red[-3], linewidth=0,whis=0,showfliers=False)
    
    # Tweak the visual presentation
    ax.xaxis.grid(True)
    ax.set(ylabel="")
    sns.despine(trim=True, left=True)
    f.savefig(figs_output+"enrichment_sig.pdf", format='pdf')
    f.savefig(figs_output+"enrichment_sig.png", format='png')   

def plot_HeritabilityPrevalence(group,tab_significant,tab_prevalence_file, figs_output, tables_output):
    
    """
    Relationships between prevalence, observed heritability estimate and number of significant genes
    Plots saved in the selected folder.
    """
    
    with open(tab_prevalence_file,"r") as f:
        tab_prevalence=pd.read_table(f,sep=",")
        
    tab_prevalence=tab_prevalence.merge(tab_significant, on="Malignancy").copy()
    tab_prevalence["no_significant_genes_99"]=[float(a) for a in tab_prevalence.no_significant_genes_99.values]
    
    nomin=tab_prevalence["prevalence"]*(1-tab_prevalence["prevalence"])
    denom=scipy.stats.norm.pdf( (scipy.stats.norm.ppf(1-tab_prevalence["prevalence"]))  )**2
    tab_prevalence["liability"]=tab_prevalence["mi_median"]*nomin/denom
    
    tab_prevalence=tab_prevalence.drop("Unnamed: 0",axis=1)
    tab_prevalence.to_csv(tables_output+"tab_liability.csv", index=False)

    plot_significant(group, tables_output, figs_output, tab_prevalence)

    fig,axes=plt.subplots(1)
    
    rp = sns.regplot(x="prevalence",y="mi_median",scatter=False,truncate=True,
                       data=tab_prevalence,ax=axes,color=palette_binary[1])
    sc=sns.scatterplot(x="prevalence",y="mi_median",data=tab_prevalence,ax=axes,**{"c":palette_binary[1]})#,size="no_significant_genes_99"

    axes.set_ylabel("median mi (h2 observed)")
    axes.set_xlabel("prevalence")

    # add annotations one by one with a loop
    for line in tab_prevalence.sort_values(by="mi_median",ascending=False).iloc[0:6].index:
        #print(tab_prevalence.Malignancy[line])
        #axes.text(tab_prevalence.prevalence[line]+0.2, tab_prevalence.mi_median[line], tab_prevalence.Malignancy[line], horizontalalignment='left', size='medium', color='black', weight='semibold')
        axes.annotate(tab_prevalence.Malignancy[line], xy=(tab_prevalence.prevalence[line], tab_prevalence.mi_median[line]),xytext=(tab_prevalence.prevalence[line]+0.0005, tab_prevalence.mi_median[line]-0.0001),rotation=-20)
    
    fig.savefig(figs_output+"prevalence.pdf", format='pdf')
    fig.savefig(figs_output+"prevalence.png", format='png')

    fig,axes=plt.subplots(1)
    
    rp = sns.regplot(x="mi_median",y="no_significant_genes_99",scatter=False,truncate=True,
                       data=tab_prevalence,ax=axes,color=palette_binary[1],ci=95)

    
    sc=sns.scatterplot(x="mi_median",y="no_significant_genes_99",sizes=(0, 1),
                       size="prevalence",data=tab_prevalence,ax=axes,**{"c":color_blu})

    axes.set_ylabel("number of significant genic regions")
    axes.set_xlabel("median mi (h2 observed)")

    # add annotations one by one with a loop
    for line in tab_prevalence.sort_values(by="mi_median",ascending=False).iloc[0:6].index:
        #print(tab_prevalence.Malignancy[line])
        #axes.text(tab_prevalence.prevalence[line]+0.2, tab_prevalence.mi_median[line], tab_prevalence.Malignancy[line], horizontalalignment='left', size='medium', color='black', weight='semibold')
        axes.annotate(tab_prevalence.Malignancy[line], xy=(tab_prevalence.mi_median[line], tab_prevalence.no_significant_genes_99[line]),xytext=(tab_prevalence.mi_median[line]+0.0005, tab_prevalence.no_significant_genes_99[line]-0.0001),rotation=-20)

    fig.savefig(figs_output+"mi_vs_noGenes.pdf", format='pdf')
    fig.savefig(figs_output+"mi_vs_noGenes.png", format='png')  

def plot_minSNP(figs_output,minSNP_file):
    
    """
    Number of minSNPs genes, HCG and overlaps between the two per cancer.
    Plots saved in the selected folder.
    """
    
    with open(minSNP_file,"r") as f:
        table=pd.read_table(f,sep=",")

    print(table.head())  

    table['tot_sigSNP']=table["sigSNP"]+table["sigSNP_HG"]
    table=table.sort_values(by=["tot_sigSNP", 'tot'], ascending=False)

    HG = table["HG"].values.tolist() 
    sigSNP = table["sigSNP"].values.tolist()
    intersection = table["sigSNP_HG"].values.tolist()

    ind = np.arange(len(table))    # the x locations for the groups

    fig,axes=plt.subplots(1, figsize=(10,8))
    fig.subplots_adjust(bottom=0.4, top=0.99)
    
    p1 = plt.bar(ind, HG,color='#0571b0')
    p3 = plt.bar(ind, sigSNP,
                 bottom=np.array(HG)+np.array(intersection), color='#ca0020')
    p2 = plt.bar(ind, intersection,
                 bottom=np.array(HG),color='#f4a582')
    
    plt.ylabel('Number of genic regions', fontsize = 8)
    #plt.title('Significant SNPs & HG')
    plt.xticks(ind, table.Malignancy.str.replace('_',' ').values.tolist(),rotation=90, fontsize = 8)
    
    plt.legend((p1[0], p2[0], p3[0]), ('HGR', 'HGR & minSNP',"minSNP"))
    

    plt.savefig(figs_output+"fig_minSNP.pdf", format='pdf')
    plt.savefig(figs_output+"fig_minSNP.png", format='png')
    
    table_normalised=table.copy()
    table_normalised["HG"]=table_normalised["HG"]/table["tot"]
    table_normalised["sigSNP"]=table_normalised["sigSNP"]/table["tot"]
    table_normalised["sigSNP_HG"]=table_normalised["sigSNP_HG"]/table["tot"]
    
    HG = table_normalised["HG"].values.tolist() 
    sigSNP = table_normalised["sigSNP"].values.tolist()
    intersection = table_normalised["sigSNP_HG"].values.tolist()
    
    ind = np.arange(len(table))    # the x locations for the groups
    width = 0.35       # the width of the bars: can also be len(x) sequence

    fig,axes=plt.subplots(1)
    p1 = plt.bar(ind, HG,color='#2c7bb6')
    p3 = plt.bar(ind, sigSNP,
                 bottom=np.array(HG)+np.array(intersection), color='#abd9e9')
    p2 = plt.bar(ind, intersection,
                 bottom=np.array(HG),color='#d7191c')
    
    plt.ylabel('percentage')
    #plt.title('Significant SNPs & HG')
    plt.xticks(ind, table.Malignancy.str.replace('_',' '),rotation=90)
    
    plt.legend((p1[0], p2[0], p3[0]), ('HG_only', 'intersection',"minSNP_only"))
    
    plt.savefig(figs_output+"fig_minSNP_normalised.pdf", format='pdf')
    plt.savefig(figs_output+"fig_minSNP_normalised.png", format='png')

    ##### TODO: breakdown
    ### breakdown
    table_breakdown=table[table['tot_sigSNP']>0]
    table_breakdown=table_breakdown.sort_values(by=['tot_sigSNP', 'tot'], ascending=True)

    NC=[int("NonCoding" in row["NotBaghera"]) for i,row in table_breakdown.iterrows()]
    table_breakdown['NC']=NC

    sigSNP = (table_breakdown["sigSNP"]-table_breakdown['NC'])#/table_breakdown['tot_sigSNP']
    NonCoding = table_breakdown['NC']#/table_breakdown['tot_sigSNP']
    intersection = table_breakdown["sigSNP_HG"]#/table_breakdown['tot_sigSNP']

    ind = np.arange(len(table_breakdown))    # the x locations for the groups

    fig,axes=plt.subplots(1)
    fig.subplots_adjust(left=0.3, right=0.99)
    
    p1 = axes.barh(ind, NonCoding,color=palette_sequential_red[-3])
    p3 = axes.barh(ind,  sigSNP,
                 left=np.array(NonCoding)+np.array(intersection), color=palette_sequential_red[7])
    p2 = axes.barh(ind, intersection,
                 left=np.array(NonCoding),color=palette_sequential_red[13])
    
    plt.xlabel('Number of Genics regions')
    #plt.title('Significant SNPs & HG')
    plt.yticks(ind, table_breakdown.Malignancy.str.replace('_',' ').values.tolist(), fontsize = 8)
    
    plt.legend((p1[0], p2[0], p3[0]), ('non-coding minSNP', 'minSNP & HGR',"minSNP"))
    
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['bottom'].set_visible(False)

    plt.savefig(figs_output+"fig_minSNP_breakdown.pdf", format='pdf')
    plt.savefig(figs_output+"fig_minSNP_breakdown.png", format='png')
    
def plot_significant(group,tables_output, figs_output, table):


    tot_SNPs=np.sum(group[group["Malignancy"]=="breast"].loc[:,["SNPs"]].values)
    
    m=[]
    summa=[]
    perc=[]
    for mal, tab in group[group['P']>0.99].groupby("Malignancy"):
        m.append(mal)
        summa.append(np.sum( tab["bg_median"]/tot_SNPs*tab["SNPs"] ))
        perc.append(100*np.sum( tab["bg_median"]/tot_SNPs*tab["SNPs"] )/tab['mi_median'].values[0])
    df_perc=pd.DataFrame()
    df_perc['Malignancy']=m
    df_perc['Sum_significant']=summa
    df_perc['Perc_significant']=perc  



    df_perc.to_csv(tables_output+'percent_significant.csv')

    table=pd.merge(table, df_perc, on='Malignancy')

    table=table.sort_values(by=['no_significant_genes_99'], ascending=False)
    table['total_perc']=100

    fig,axes=plt.subplots(1, 3, figsize=(8,10),gridspec_kw={'width_ratios': [1,3,4],"wspace":0.025})
    fig.subplots_adjust(left=0.2,right=0.80) 
    
    vec=np.asarray(table.mi_median.values)[:,np.newaxis]
    g=sns.heatmap(vec,cmap=palette_sequential_blu,ax=axes[0],xticklabels=0,yticklabels=1,fmt="0.4f",annot=True,cbar=False,annot_kws={"size": 6})
    g.set_yticklabels(table["Malignancy"].str.replace('_', " ") , rotation = 0, fontsize = 8)
    g.set_xlabel("h2 obs", rotation = 0, fontsize = 10) 

    #vec=np.asarray(table.Perc_significant.values)[:,np.newaxis]
    #g=sns.heatmap(vec,cmap=palette_sequential_blu,ax=axes[1],xticklabels=1,yticklabels=0,fmt="0.3f",annot=True,cbar=False,annot_kws={"size": 6})
    #g.set_xticklabels(["h2_liab"], rotation = 0, fontsize = 8) 

    #point=sns.pointplot(x="Perc_significant", y="Malignancy",
    #          data=table, color=palette_binary[0],join=False,
    #          markers="o", ax=axes[1])
    sns.barplot(x="total_perc", y="Malignancy", data=table,
                color='#fddbc7',ax=axes[1])
    sns.barplot(x="Perc_significant", y="Malignancy", data=table,
                color='#b2182b',ax=axes[1])



    axes[1].set_yticklabels("")
    axes[1].set( ylabel="",
           xlabel="%h2 explained by HG")
    #point.legend("")
    #axes[1].grid(b=True, which='minor', axis='x', lw=0.1, linestyle='--')
    axes[1].legend().set_visible(False)

    #sns.barplot(x="no_significant_genes_95", y="Malignancy", data=table,
    #            label="P > 0.95", color=palette_sequential_blu[4],ax=axes[1])
                
    #m=0            
    #for p in ax.patches:
        #ax.annotate( tab3.mi_median.values.tolist()[m]  , (p.get_width(), p.get_y()+0.5 ),fontsize=6,horizontalalignment='left')
        #m+=1
    
    bar=sns.barplot(x="no_significant_genes_99", y="Malignancy", data=table,
                label="%s > 0.99" %str(chr(951)), color=palette_sequential_blu[-1],ax=axes[2])
    bar.set_yticklabels("")
    # Add a legend and informative axis label
    axes[2].legend(ncol=1, loc="lower right", frameon=True)
    axes[2].set( ylabel="",
           xlabel="Number of heritability genic regions")

    axes[2].set_axisbelow(True)
    axes[2].xaxis.grid(color='gray', linestyle='dotted',linewidth=0.5)
    
    sns.despine(left=True, bottom=True)
    fig.savefig(figs_output+"fig_significant.pdf", format='pdf')
    fig.savefig(figs_output+"fig_significant.png", format='png')    



#plot_minSNP('/Users/s1899202/Documents/','/Users/s1899202/Documents/Projects/baghera_project/results/final_results/tables/SNPvsHG.csv')