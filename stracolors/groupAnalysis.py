#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 17:23:48 2018

@author: viola
"""

#Group Analysis of Cancer

import matplotlib as mpl
mpl.use("TkAgg")

import numpy as np
import pandas as pd
import gseapy as gp
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

import figure1
import figure2
import figure3


def reset_matplotlib():
    """
    Reset matplotlib to a common default.
    """
    # Set all default values.
    mpl.rcdefaults()
    # Force agg backend.
    plt.switch_backend('TkAgg')
    #matplotlib.use("TkAgg")
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    mpl.rcParams['font.family'] = 'Arial'

    #mpl.rcParams['font.family'] = 'serif'
    #mpl.rcParams['font.serif'] = ['Charter']

    mpl.rcParams['axes.spines.left'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.bottom'] = False
    mpl.rcParams['axes.spines.top'] = False

    mpl.rcParams['axes.facecolor']= 'white'
    mpl.rcParams['axes.edgecolor']= '#264258'


    mpl.rcParams['grid.alpha']=0.5
    mpl.rcParams['grid.color']='#264258' ## grid color
    mpl.rcParams['grid.linestyle']='dotted'        ## solid
    mpl.rcParams['grid.linewidth']= 0.5      ## in points

    
    


#font_manager.findSystemFonts()
#matplotlib.font_manager._rebuild()
#prop = font_manager.FontProperties(fname="/usr/share/fonts/truetype/lato/Lato-Bold.ttf")
#matplotlib.rcParams['font.family'] = prop.get_name()
#matplotlib.rcParams['font.weight']="bold"

#def parserFunc():
#    parser = argparse.ArgumentParser(description='Parse command line options.')
#    parser.add_argument('--folderInput', '-f',  type=str, required=True, help="Data Input, use the genesResults.csv file, output of regression.py ")
#    parser.add_argument('--pathwaySet', '-p',  type=inType, required=True, help="Type of pathway set KEGG or canonical")
#    return parser

def weightedTtest(table, column, genelist_aggregate=None,genelist_symbol=None):
    '''
    one-tailed ttest of the values in a column, testing genes in the genelist
    against those outside. Returns stats and pvalue
    '''
    table=table.reset_index()
    IN=0
    OUT=0
    ind=[]
    if genelist_symbol:

        for g in genelist_symbol:                      
            ind.append(table[table['name'].str.contains("^"+g+"[;]|^"+g+"$|[;] "+g+"$|[;] "+g+"[;]")].index.values)
            
        ind=list(filter(None, ind))
        print(ind)

        #TODO Change to set
        ind = list(set([int(x) for x in ind]))

    else:
        for g in genelist_aggregate:                      
            ind.append(table[table['name']==g].index.values)
        ind=list(filter(None, ind))
        #TODO Change to set
        ind = list(set([int(x) for x in ind]))
        
    IN=table[column].iloc[ind]
    OUT=table[column].iloc[list(set(table.index).difference(set(ind)) )]
            
    stat,p=stats.ttest_ind(IN,OUT)
    if stat>0:
        pvalue=p/2
    else:
        pvalue=1     

    return stat,pvalue

def hallmarksOverlaps(group, cosmic_file):
    
    """This function finds the overlapping genes between the heritability analysis significant ones and
    those in the cosmic hallmarks 
    """

    with open(cosmic_file, "r") as f:
        COSMIC=pd.read_table(f, sep=',')

    #only somatic mutations
    COSMIC=COSMIC[COSMIC["Somatic"]=="yes"]
    
    COSMIC_hallmark=COSMIC[COSMIC["Hallmark"]=="Yes"]
    COSMIC_hallmark_set=set(COSMIC_hallmark["Gene Symbol"].values.tolist())

    ALL=[]
    for gene,gr in group.groupby("name"):
        for i in gene.split(';'):
            ALL.append(i.replace(" ", ""))
            
    P_hallmarks=set(ALL).intersection(COSMIC_hallmark_set)
        
    overlaps=group[group["P"]>0.99].copy()   
    
    all_significant=[]
    for gene,gr in overlaps.groupby("name"):
        for i in gene.split(';'):
            all_significant.append(i.replace(" ", ""))
    
    P_baghera=len(set(all_significant))
    
    PP_set=set(all_significant).intersection(P_hallmarks)
    PP=len(PP_set)
    PP_list=list(PP_set)
    NP=len(P_hallmarks)-PP
    PN=P_baghera-PP
    NN=len(ALL)-PP-NP-PN
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")
    print("OR fisher-test hallmarks= "+str(oddsratio))
    print("pvalue fisher-test hallmarks= "+str(pvalue))  
    
    with open(stats_file, "a") as f:    
        f.write("\n \n---------Hallmarks Analysis----------") 

    with open(stats_file, "a") as f:          
        f.write("\n \n #### Hallmarks enrichment")   
        f.write("\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d \t" %(PP,PN,NP,NN) )        
        f.write("\n OR: %f, Fisher p_value: %f \t" %(oddsratio, pvalue))
        f.write("\n number PP: "+str(PP))
        f.write("\t\n PP:"+str(PP_list))




def DNArepairOverlaps(group, DNArepair_file):
    
    """This function finds the overlapping genes between the heritability analysis significant ones and
    those in the DNA repair dataset 
    """

    with open(DNArepair_file, "r") as f:
        repair=pd.read_table(f)
        
    n_repair=len(repair)
    
    repair_list=repair.iloc[:,0].str.replace("(","").str.replace(")","").replace(","," ").str.split(" ").values.tolist()
    repair_set=set([j[0] for j in repair_list])
    repair_set_synonyms=set([j for i in repair_list for j in i[1:] if j!=""])
    repair_set_total=set([j for i in repair_list for j in i if j!=""])

    ALL=[]
    for gene,gr in group.groupby("name"):
        for i in gene.split(';'):
            ALL.append(i.replace(" ", ""))
            
    P_repair=list(set(ALL).intersection(repair_set))    
    P_repair_s=list(set(ALL).intersection(repair_set_synonyms))    
        
    overlaps=group[group["P"]>0.99].copy()   
    
    all_significant=[]
    for gene,gr in overlaps.groupby("name"):
        for i in gene.split(';'):
            all_significant.append(i.replace(" ", ""))
    
    P_baghera=len(set(all_significant))
    
    PP_set=set(all_significant).intersection(repair_set_total)
    PP=len(PP_set)
    PP_list=list(PP_set)
    NP=len(P_repair)-PP
    PN=P_baghera-PP
    NN=len(ALL)-PP-NP-PN
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")
    print("OR fisher-test DNA repair= "+str(oddsratio))
    print("pvalue fisher-test DNA repair= "+str(pvalue))  

    oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False)       
    #oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")        
    
    print("OR chi2 DNA repair= "+str(oddsratio))
    print("pvalue chi2 DNA repair= "+str(pvalue))      
    
    with open(stats_file, "a") as f:    
        f.write("\n \n---------DNA repair Analysis----------") 

    with open(stats_file, "a") as f:          
        f.write("\n \n #### DNA repair enrichment")   
        f.write("\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d \t" %(PP,PN,NP,NN) )        
        f.write("\n OR: %f, Fisher p_value: %f \t" %(oddsratio, pvalue))
        f.write("\n number PP: "+str(PP))
        f.write("\t\n PP:"+str(PP_set))




def GOOverlaps(group):
    
    """This function find the overlapping genes between the heritability analysis significant ones and
    those in the COSMIC dataset reported as genes with driver somatic mutations 
    """
    #creates the gmt file of the overlaps between BAGHERA and COSMIC
#    gmt_dict={}
#    gmt_dict["somatic_overlaps"]={}
#    gmt_dict["non_overlapping"]={}
#    gmt_dict["somatic_overlaps"]["descriptor"]="COSMIC_BAGHERA"    
#    gmt_dict["non_overlapping"]["descriptor"]="BAGHERA_no_COSMIC"
#    gmt_dict["all_significant"]={}
#    gmt_dict["all_significant"]["descriptor"]="ALL_BAGHERA"    
    # All significant genes     

    ALL=[]
    for gene,gr in group.groupby("name"):
        for i in gene.split(';'):
            ALL.append(i.replace(" ", ""))
            
    GOtot=set(GO["gene"].values )
    total=list(set(ALL).intersection(GOtot))     
    len(total)

        
    overlaps=group[group["P"]>0.99].copy()   
    
    all_significant=[]
    for gene,gr in overlaps.groupby("name"):
        for i in gene.split(';'):
            all_significant.append(i.replace(" ", ""))
    
    n_sig=len(all_significant)
    
    all_significant=set(all_significant).intersection(total)
    
    with open(stats_file, "a") as f:    
        f.write("\n \n---------Gene Ontology Analysis----------")                   
        f.write("\n \nThere are %d genic areas in total " %len(set(group["name"].values)) )            
        f.write("\n \nThere are %d genes in total " %len(ALL))   
        f.write("\n \nThere are %d genes in GO in total" %len(GOtot))   
        f.write("\n \nOf those, %d genes are in our geneset " %len(total))
        f.write("\n \nThere are %d significant genes " %n_sig)
        f.write("\n \nOf those, %d genes are in GO" %len(all_significant))
    
    GOintersect=GO[GO["gene"].isin(total)]    
    
    goid=[]
    OR=[]
    p=[]
    df=pd.DataFrame()
    PP_names=[]
    for GO_id,tab in GOintersect.groupby(["goid"]):
        goid.append(GO_id)
        print(GO_id)
        print(len(tab))
                
        PP=len(set(tab["gene"]).intersection(all_significant))
        PP_list=list(set(tab["gene"]).intersection(all_significant))
        NP=len(set(tab["gene"]).intersection(total))-PP
        PN=len(all_significant)-PP
        NN=len(total)-PP-NP-PN
        oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
        #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
        OR.append(oddsratio)
        p.append(pvalue)
        PP_names.append(PP_list)
        
             
    
    df["id"]=goid
    df["OR"]=OR
    df["pvalue"]=p   
    df["PP_list"]=PP_names
    
    rejects,pval,k,bonf=smm.multipletests(df['pvalue'],alpha=0.1,method="fdr_bh")
    df["rejects"]=rejects
    df["bh_pvalue"]=pval
    df["k"]=k
    df["bonf"]=bonf
    
    
    GOdrop=GOintersect[["goid","term","ontology"]].drop_duplicates()

    df=df.merge(GOdrop,how="inner", left_on='id', right_on='goid')
    df.to_csv(tables_output+"GO_overlaps.csv")
    
    
    return df

def KEGGOverlaps(group):
    
    """This function find the overlapping genes between the heritability analysis significant ones and
    those in the COSMIC dataset reported as genes with driver somatic mutations 
    """
    #creates the gmt file of the overlaps between BAGHERA and COSMIC
#    gmt_dict={}
#    gmt_dict["somatic_overlaps"]={}
#    gmt_dict["non_overlapping"]={}
#    gmt_dict["somatic_overlaps"]["descriptor"]="COSMIC_BAGHERA"    
#    gmt_dict["non_overlapping"]["descriptor"]="BAGHERA_no_COSMIC"
#    gmt_dict["all_significant"]={}
#    gmt_dict["all_significant"]["descriptor"]="ALL_BAGHERA"    
    # All significant genes
    
    ALL=[]
    for gene,gr in group.groupby("name"):
        for i in gene.split(';'):
            ALL.append(i.replace(" ", ""))
    
    KEGGtot=[]
    for k,val in KEGG.items():
        KEGGtot=list(set(KEGGtot).union(set(val)))
    
    total=list(set(ALL).intersection(KEGGtot))     
    len(total)
        
    overlaps=group[group["P"]>0.99].copy()   
    
    all_significant=[]
    for gene,gr in overlaps.groupby("name"):
        for i in gene.split(';'):
            all_significant.append(i.replace(" ", ""))
    
    all_significant=set(all_significant).intersection(total)
    
    #GOintersect=GO[GO["gene"].isin(total)]    
    
    paths=[]
    OR=[]
    p=[]
    df=pd.DataFrame()
    for path,vals in KEGG.items():
        paths.append(path)
        #print(path)
        IN=set(vals).intersection(total)
                
        PP=len(IN.intersection(all_significant))
        NP=len(IN)-PP
        PN=len(all_significant)-PP
        NN=len(total)-PP-NP-PN
        #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
        oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
        OR.append(oddsratio)
        p.append(pvalue)
        
        #print("pvalue fisher-test= "+str(pvalue))      
    
    df["path"]=paths
    df["OR"]=OR
    df["pvalue"]=p    
    
    rejects,pval,k,bonf=smm.multipletests(df['pvalue'],alpha=0.1,method="fdr_bh")
    df["rejects"]=rejects
    df["bh_pvalue"]=pval
    df["k"]=k
    df["bonf"]=bonf
 
    df.to_csv(tables_output+"KEGG_overlaps.csv")   
    return df

def createTablePivot(group, genes_of_interest, setname):
    
    """
    Computes the heritability explained by a list of genes. 
    It returns a dictionary with three pivot tables, 
        - one for all genes
        - one for significant genes
        - one for panCancer genes
    """

    GOI_set=set(genes_of_interest)
    with open(stats_file, "a") as f:          
        f.write("\t\n \n Pivot Table for %s genes\n-------------------------- " %(setname))   
        f.write("\t\n \n there are %d genes in the set" %(len(GOI_set)))   

    # Number of SNPs
    tot_SNPs=np.sum(group[group["Malignancy"]=="breast"].loc[:,["SNPs"]].values)


    group["bg_weighted"]=group["SNPs"]*group["bg_median"]/tot_SNPs
    
    mi={}
    mi["bg_tot"]={}
    mi["mi"]={}
    for malignancy,tab in group.groupby(["Malignancy"]):
        mi["bg_tot"][malignancy]=np.sum( tab["bg_median"]/tot_SNPs*tab["SNPs"] )
        mi["mi"][malignancy]=tab.iloc[0].loc["mi_median"]
    
    
    pivot={}

    for th in [0,1,2]:

        if th==0:
            sig_genes=[name for name,g in group.groupby(["name"])]
        else:
        # Check which genes are significant in more than 0,1,2 malignancies
            sig_genes=[name for name,g in group[group["P"]>0.99].groupby(["name"]) if len(g)>=th]

        overlaps=[]
        for k in GOI_set:
            for j in sig_genes:
                if k in j.replace(" ","").split(';'):
                    overlaps.append(j)
        
        if len(overlaps)>1:

            with open(stats_file, "a") as f:   
                f.write("\t\n \n **Threshold>=%d** ttest with pval<0.1:" %(th) )          
        
            tab_ttest=pd.DataFrame()
            for mal, tab in group.groupby('Malignancy'):
                s, pval= weightedTtest(tab, 'bg_weighted', genelist_aggregate=overlaps)
                if pval<0.1:
                    with open(stats_file, "a") as f:   
                        f.write("\t\nFor %s : stats= %f, p-value= %f" %(mal, s, pval))   
                tab_ttest=tab_ttest.append({'Malignancy': mal,
                                        'len_genelist': len(overlaps), 
                                        'stats':s ,
                                        'pval':pval}, ignore_index=True)

            tab_ttest.to_csv(tables_output+"ttest_"+setname+"th>="+str(th)+".csv", sep=",", mode='w',index=True)
        

            with open(stats_file, "a") as f:   
                f.write("\t\n \n **Threshold>=%d** (in the table the genes are significant in at least %d malignancies)" %(th,th) )          
                f.write("\t\n \n There are %d genes in the table" %(len(sig_genes) ) )   
                f.write("\t\n \n There are %d genes in the table that are %s" %(len(overlaps),setname ) )   

            # Table with the genes that are both significant and in the GOI set
            tableGOI=group[group["name"].isin(overlaps)]

            # Create the pivot table
            pivotGOI = tableGOI.pivot("name", "Malignancy", "bg_weighted")
            n_genes=len(pivotGOI)

            # Count the number of SNPs in the GOI
            SNPs=[]
            for i,tab in pivotGOI.iterrows():
                SNPs.append(tableGOI[tableGOI["name"]==i]["SNPs"].values[0])
            SNPs.append(tot_SNPs)
            SNPs.append(sum(SNPs[:-1]))
            SNPs.append(tot_SNPs)
            SNPs.append(SNPs[-2]/tot_SNPs)
            
            mal=pivotGOI.columns
            
            # Add rows for mi_median and bg_tot to the pivot table
            mi_list=[tableGOI[tableGOI["Malignancy"]==m].iloc[0].loc["mi_median"] for m in mal]    
            bg_list=[mi["bg_tot"][m] for m in mal] 
            df2 = pd.DataFrame([np.array(mi_list),np.array(bg_list)], columns=pivotGOI.columns)
            #df2["ind"]="mi_median"
            #df2.set_index(["ind"])
            df2=df2.rename({0: "mi_median", 1:"bg_tot"})
            pivotGOI= pivotGOI.append(df2, ignore_index=False)
            
            #pivot_panGenes["av_weighted"]=np.average(pivot_panGenes.values,axis=1)
            
            sum_mal=np.sum(pivotGOI.values[:-2,:],axis=0)
            total_ratio=sum_mal/np.array(mi_list)
            
            df3 = pd.DataFrame([sum_mal.T,total_ratio.T,], columns=pivotGOI.columns)
            df3=df3.rename({0: "mi_median",1: "ratio"})
            pivotGOI= pivotGOI.append(df3, ignore_index=False)

            # The column needs to be added now, but was created before
            pivotGOI["SNPs"]=SNPs
                
            pivot[th]=pivotGOI

            pivotGOI.to_csv(tables_output+"pivot"+setname+"th>="+str(th)+".csv", sep=",", mode='w',index=True)
        else:
            with open(stats_file, "a") as f:   
                f.write("\t\n \n **Threshold>=%d**, there are less than 2 overlaps , analysis skipped " %(th) )   
                f.write("\t\n Overlaps:"+str(overlaps))   
            pivot[th]=pd.DataFrame()

    return pivot

def panCancer_binary(group):
    
    """
    Creates pivot tables for panCancers
    """
    tot_SNPs=np.sum(group[group["Malignancy"]=="breast"].loc[:,["SNPs"]].values)
    
    significant_list_splitted=list(set([i for j in group[group["P"]>0.99].name.str.replace(" ","").str.split(";") for i in j]))
    
    pivot={}
    for th in [0,1]:
        
        sig_genes=[name for name,g in group[group["P"]>0.99].groupby(["name"]) if len(g)>th]
        #sig_genes=list(set(group[group["P"]>0.99]["name"].values))   
           
        panGenes=group[group["name"].isin(sig_genes)]
        
        panGenes['binary']=1.*panGenes['P']>0.99

        print(panGenes.head())
        pivot_panGenes = panGenes.pivot("name", "Malignancy", "binary")
            
        pivot[th]=pivot_panGenes 

        pivot_panGenes.to_csv(tables_output+"binary_panGenes_th>"+str(th)+".csv", sep=",", mode='w',index=True)


    return pivot




def panCancerEnrichment(group):
    
    """
    Creates pivot tables for panCancers
    """
    tot_SNPs=np.sum(group[group["Malignancy"]=="breast"].loc[:,["SNPs"]].values)
    
    
    mi={}
    mi["bg_tot"]={}
    mi["mi"]={}
    for malignancy,tab in group.groupby(["Malignancy"]):
        print(malignancy)
        mi["bg_tot"][malignancy]=np.sum( tab["bg_median"]/tot_SNPs*tab["SNPs"] )
        mi["mi"][malignancy]=tab.iloc[0].loc["mi_median"]
    
    significant_list_splitted=list(set([i for j in group[group["P"]>0.99].name.str.replace(" ","").str.split(";") for i in j]))
    
    pivot={}
    for th in [0,1]:
        
        sig_genes=[name for name,g in group[group["P"]>0.99].groupby(["name"]) if len(g)>th]
        #sig_genes=list(set(group[group["P"]>0.99]["name"].values))   
           
        panGenes=group[group["name"].isin(sig_genes)]
        
        panGenes["bg_weighted"]=panGenes["SNPs"]*panGenes["bg_median"]/tot_SNPs
        panGenes['binary']=panGenes['P']>0.99

        print(panGenes.head())
        pivot_panGenes = panGenes.pivot("name", "Malignancy", "bg_weighted")

        
        n_genes=len(pivot_panGenes)
        
        mal=pivot_panGenes.columns
        
        SNPs=[]
        for i,tab in pivot_panGenes.iterrows():
            SNPs.append(panGenes[panGenes["name"]==i]["SNPs"].values[0])
        SNPs.append(tot_SNPs)
        SNPs.append(sum(SNPs[:-1]))
        SNPs.append(tot_SNPs)
        SNPs.append(SNPs[-2]/tot_SNPs)
        
               
        mi_list=[panGenes[panGenes["Malignancy"]==m].iloc[0].loc["mi_median"] for m in mal]    
        bg_list=[mi["bg_tot"][m] for m in mal] 
        df2 = pd.DataFrame([np.array(mi_list),np.array(bg_list)], columns=pivot_panGenes.columns)
        #df2["ind"]="mi_median"
        #df2.set_index(["ind"])
        df2=df2.rename({0: "mi_median", 1:"bg_tot"})
        pivot_panGenes= pivot_panGenes.append(df2, ignore_index=False)
        
        #pivot_panGenes["av_weighted"]=np.average(pivot_panGenes.values,axis=1)
        
        sum_mal=np.sum(pivot_panGenes.values[:-2,:],axis=0)
        total_ratio=sum_mal/np.array(mi_list)
        
        df3 = pd.DataFrame([sum_mal.T,total_ratio.T,], columns=pivot_panGenes.columns)
        df3=df3.rename({0: "mi_median",1: "ratio"})
        pivot_panGenes= pivot_panGenes.append(df3, ignore_index=False)
        

        panGenes_list_splitted=list(set([i for j in sig_genes for i in j.replace(" ","").split(";")]))
        
            
        print("%d genes are significant in at least %d cancers" %(len(panGenes_list_splitted),th+1))
        
        pivot_panGenes["SNPs"]=SNPs
            
        pivot[th]=pivot_panGenes 
        
        #pivot_panGenes= pivot_panGenes.append(pivot_panGenes.sum(numeric_only=True)/n_genes, ignore_index=True)
        #pivot_panGenes= pivot_panGenes.append(pivot_panGenes.sum(numeric_only=True)/pivot_panGenes.iloc[-2], ignore_index=False)
    
        
        pivot_panGenes.to_csv(tables_output+"panGenes_th>"+str(th)+".csv", sep=",", mode='w',index=True)

        pivot_binary_panGenes = panGenes.pivot("name", "Malignancy", "binary")
    
    return pivot


     
        
def panCancerTable(group):
    
    """
    Creates pan cancer table
    """
    sig_table=group[group["P"]>0.99].copy()
       
    panCancer=pd.DataFrame()
    
    for genes, tab in sig_table.groupby(["name"]):
        
        n_cancer_ukbb=len(tab["Malignancy"])
        name_cancer_ukbb=tab["Malignancy"].values.tolist()
        min_enrichment=np.min(tab["weight"].values)
        max_enrichment=np.max(tab["weight"].values)
        avg_enrichment=np.average(tab["weight"].values)
        
        
        
        
        panCancer=panCancer.append({'name': genes,
                                    'n_cancer_ukbb': n_cancer_ukbb, 
                                    'name_cancer_ukbb':name_cancer_ukbb,
                                    'min_enrichment':min_enrichment,
                                    'max_enrichment':max_enrichment,
                                    'avg_enrichment':avg_enrichment}, ignore_index=True)
    
    panCancer=annotateGeneTable(panCancer.copy())
    panCancer=annotateKEGG(panCancer.copy())
    
    panCancer=panCancer.sort_values(by=["n_cancer_ukbb"],ascending=False)
    
    panCancer.to_csv(tables_output+"panCancerTable.csv", sep=",", mode='w',index=False)
    
    return panCancer
    
def table_significant(group):
    """Figure with no_significant genes vs no of cancers 
    """
    
    tab3=pd.DataFrame(columns=['Malignancy','no_significant_genes_95','no_significant_genes_99', 'Top5(P_and_bg_median )'])    
    tab3_complete=pd.DataFrame(columns=['Malignancy','no_significant_genes_95','no_significant_genes_99', 'Top5(P_and_bg_median )'])    
    
    #signif=group[group["Significant"]==True].copy()
    for mal,g in group.groupby('Malignancy'):
        g=g.sort_values(by=["P","bg_median"])
        mi=g.mi_median.iloc[0]
        tab3 = tab3.append({'Malignancy': mal,'mi_median':mi,'no_significant_genes_95':len(g[g["P"]>0.95]),'no_significant_genes_99':len(g[g["P"]>0.99]),  'Top5(P_and_bg_median )':g.iloc[-5:].name.values}, ignore_index=True)
        tab3_complete = tab3_complete.append({'Malignancy': mal,'mi_median':mi,'mi_5':mi,'mi_95':mi,'no_significant_genes_95':len(g[g["P"]>0.95]),'no_significant_genes_99':len(g[g["P"]>0.99]),  'Top5(P_and_bg_median )':g.iloc[-5:].name.values}, ignore_index=True)
                  
    tab3.to_csv(tables_output+"tab3_signif.csv", sep=",", mode='w')
    
    tab3=tab3.sort_values(by="no_significant_genes_95",ascending=False)
    
    return tab3

    
def actionableOverlaps(group):
    
    """This function find the overlapping genes between the heritability analysis significant ones and
    those in the cancerGenesList
    """
    
    
    actionable_genes_list=list(set(actionable.Gene.values.tolist()))
    
    sig=group[group["P"]>0.99].copy()    
    sig["actionable"]="no"
    
    for index, row in sig.iterrows():
        for i in row["name"].split(';'):
            gene_name=str.upper(i.replace(" ", ""))
            if gene_name in actionable_genes_list:
                print(gene_name)                
                sig.set_value(index,"actionable","yes")

    sig=sig[sig["actionable"]=="yes"]
    tab_actionable = sig.pivot("name", "Malignancy", "bg_median")


    fig,axes=plt.subplots(1, 1, figsize=(10, 10),sharey=True)
    fig.subplots_adjust(left=0.2,right=0.99) 
    #pal=sns.cubehelix_palette(100)
    
    # Draw a heatmap with the numeric values in each cell
    g=sns.heatmap(tab_actionable,xticklabels=1,yticklabels=1,
                  cmap=Teal_7.mpl_colormap,ax=axes,linewidths=.5, cbar_kws={"orientation": "vertical","aspect":40})
    
    g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 6)
    g.set_xticklabels(g.get_xticklabels(), rotation = 90, fontsize = 6) 
    
    g.set_title("Actionable targets")
    fig.savefig(figs_output+"actionable.pdf", format='pdf')
    fig.savefig(figs_output+"actionable.jpeg", format='jpeg')
    tab_actionable.to_csv(tables_output+"actionable.csv", sep=",", mode='w')
    


    #g=sns.clustermap(common_genes_bin.iloc[:,0:-1],method="single",cmap=pal,row_cluster=False, **{'xticklabels':1})
    #g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation = 0, fontsize = 6)
    #g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation = 90, fontsize = 6
def printGMT(dataDict, fileOutput):
    
    """function to write GMT files from dictionaries, data dict keys are the gene-sets names, for each gene set there 
    must be a dictionary with the key "descriptor" as the description name for the set and a "genes" key with the list of genes
    the output file will be located into the fileOutput path    """
        
    output=fileOutput
    with open(output, 'w') as f:
        f.write('')    
    
    for key, dict_set in dataDict.items():
        with open(output, 'a') as f:        
            f.write(str(key)+'\t'+str(dict_set["descriptor"])+'\t'+'\t'.join(dict_set["genes"])+"\n")

def annotateGeneTable(genes_table):
    """
    Annotates a table of genes (requires name as column) in the results, annotations in cancerGenesList
    """
    
    #genes_table["Significant"]=genes_table["P"]>0.99
    #genes_table["Fold_enrichment"]=(genes_table["bg_median"]-genes_table["mi_median"])/genes_table["mi_median"]
    oncokb_list=cancerGenesList[cancerGenesList["OncoKB Annotated"]=="Yes"]["Hugo Symbol"].values.tolist()
    #print("There are %d genes annotated by oncoKB" %len(oncokb_list))    
    oncogene_list=cancerGenesList[cancerGenesList["OncoKB Oncogene"]=="Yes"]["Hugo Symbol"].values.tolist()
    #print("There are %d genes annotated by oncoKB as oncogene" %len(oncogene_list))
    TSG_list=cancerGenesList[cancerGenesList["OncoKB TSG"]=="Yes"]["Hugo Symbol"].values.tolist()
    #print("There are %d genes annotated by oncoKB as TSG" %len(TSG_list))
    CGC_list=cancerGenesList[cancerGenesList["Sanger CGC"]=="Yes"]["Hugo Symbol"].values.tolist()
    #print("There are %d genes annotated by Sanger CGC" %len(CGC_list))
    all_cancer_list=cancerGenesList["Hugo Symbol"].values.tolist()
    #print("There are %d genes annotated as relevant from the oncoKB project" %len(all_cancer_list))
    
    genes_table["OncoKB_list"]=[list(set(a).intersection(set(oncokb_list))) for a in genes_table.name.str.replace(" ","").str.split(";")]
    genes_table["OncoKB"]=genes_table.apply(lambda row: len(row["OncoKB_list"]) >0, axis=1)
    
    genes_table["OncoKB_Oncogene_list"]=[list(set(a).intersection(set(oncogene_list))) for a in genes_table.name.str.replace(" ","").str.split(";")]
    genes_table["OncoKB_Oncogene"]=genes_table.apply(lambda row: len(row["OncoKB_Oncogene_list"]) >0, axis=1)
    
    genes_table["OncoKB_TSG_list"]=[list(set(a).intersection(set(TSG_list))) for a in genes_table.name.str.replace(" ","").str.split(";")]
    genes_table["OncoKB_TSG"]=genes_table.apply(lambda row: len(row["OncoKB_TSG_list"]) >0, axis=1)
    
    genes_table["CGC_list"]=[list(set(a).intersection(set(CGC_list))) for a in genes_table.name.str.replace(" ","").str.split(";")]
    genes_table["CGC"]=genes_table.apply(lambda row: len(row["CGC_list"]) >0, axis=1)
    
    genes_table["all_cancer_list"]=[list(set(a).intersection(set(all_cancer_list))) for a in genes_table.name.str.replace(" ","").str.split(";")]
    genes_table["all_cancer"]=genes_table.apply(lambda row: len(row["all_cancer_list"]) >0, axis=1)
    
    genes_table["Actionable_list"]=[list(set(a).intersection(set(actionable["Gene"].values.tolist()))) for a in genes_table.name.str.replace(" ","").str.split(";")]
    genes_table["Actionable"]=genes_table.apply(lambda row: len(row["Actionable_list"]) >0, axis=1)

    
    #genes_table=genes_table.sort_values(by=["P"],ascending=False)
    
    return genes_table

def annotateKEGG(genes_table):
    """
    Annotates a table of genes (requires name as column) in the results, annotations in KEGG
    """

    L=["" for i in range(len(genes_table))]
    
    for path, genes in KEGG.items():
        path_genes=[list(set(a).intersection(set(genes))) for a in genes_table.name.str.replace(" ","").str.split(";")]
        temp=[path+"," if len(i)>0 else "" for i in path_genes for i in path_genes]
        for i in range(len(L)):
            L[i]=L[i]+temp[i]
        
    genes_table["pathways"]=L           

    
    #genes_table=genes_table.sort_values(by=["P"],ascending=False)
    
    return genes_table

#%%          
###########################################################
###### SECTIONS ###########################################
###########################################################


def oncoKB_enrichments(group, annotated, panCHG_table):
    
    significant=group[group["P"]>0.99]
    
    all_genes_single=set([b for a in group.name.str.replace(" ","").str.split(";") for b in a])
    with open(stats_file, "a") as f:    
        f.write("\n \n ### OncoKB Enrichments" )     
        f.write("\n \nThere are %d genic areas in total \t" %len(set(group["name"].values)) )            
        f.write("\n \nThere are %d genes in total \t" %len(all_genes_single))    
    
    
    with open(stats_file, "a") as f:
        
        oncokb_list=cancerGenesList[cancerGenesList["OncoKB Annotated"]=="Yes"]["Hugo Symbol"].values.tolist()
        f.write("\n \nThere are %d genes annotated by oncoKB" %len(oncokb_list))    
        f.write("\n \nOut of the %d genes annotated by oncoKB, %d are in our geneset too" %(len(oncokb_list),len(all_genes_single.intersection(set(oncokb_list)))))    
        
        oncogene_list=cancerGenesList[cancerGenesList["OncoKB Oncogene"]=="Yes"]["Hugo Symbol"].values.tolist()
        f.write("\n \nThere are %d genes annotated by oncoKB as oncogene" %len(oncogene_list))
        f.write("\n \nOut of the %d genes annotated by oncoKB as oncogene, %d are in our geneset too" %(len(oncogene_list),len(all_genes_single.intersection(set(oncogene_list))))  )  

        TSG_list=cancerGenesList[cancerGenesList["OncoKB TSG"]=="Yes"]["Hugo Symbol"].values.tolist()
        f.write("\n \nThere are %d genes annotated by oncoKB as TSG" %len(TSG_list))
        f.write("\n \nOut of the %d genes annotated by oncoKB as TSG, %d are in our geneset too" %(len(TSG_list),len(all_genes_single.intersection(set(TSG_list))))    )

        f.write("\n \nThere are %d genes annotated as both TSG and oncogene" %len(set(TSG_list).intersection(set(oncogene_list))))
        f.write("\n \nThere are %d genes that are in oncoKB but with unknown function" %len(set(oncokb_list).difference(set(oncogene_list)).difference(set(TSG_list)))    )


        CGC_list=cancerGenesList[cancerGenesList["Sanger CGC"]=="Yes"]["Hugo Symbol"].values.tolist()
        f.write("\n \nhere are %d genes annotated by Sanger CGC" %len(CGC_list))
        f.write("\n \nOut of the %d genes annotated by Sanger CGC, %d are in our geneset too" %(len(CGC_list),len(all_genes_single.intersection(set(CGC_list))))    )

        all_cancer_list=cancerGenesList["Hugo Symbol"].values.tolist()
        f.write("\n \nThere are %d genes annotated as relevant from the oncoKB project" %len(all_cancer_list))
        f.write("\n \nOut of the %d genes annotated as relevant from the oncoKB project, %d are in our geneset too" %(len(all_cancer_list),len(all_genes_single.intersection(set(all_cancer_list)))) )   




        
    ################# Enrichment for CHG ######################################
    CHG_single=set([b for a in significant.name.str.replace(" ","").str.split(";") for b in a])
    gmt_dict["CHGs"]={}
    gmt_dict["CHGs"]["descriptor"]="baghera"
    gmt_dict["CHGs"]["genes"]=list(CHG_single)
    
    
    with open(stats_file, "a") as f:          
        f.write("\n \n #### Enrichment CHG \t")         
        
        f.write("\n \nThere are %d significant genic areas in total \t" %len(set(significant["name"].values)) )    
        f.write("\n \nThey corresponds are %d genes in total \t " %len(CHG_single) )    


    # oncoKB annotated
    tot=len(all_genes_single)   
    PP=len(CHG_single.intersection(set(oncokb_list)))
    PN=len(all_genes_single.intersection(set(oncokb_list))  )-PP
    NP=len(CHG_single)-PP
    NN=tot-PP-PN-NP
    oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    #oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\n \n #### OncoKB annotated enrichment")   
        f.write("\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d \t" %(PP,PN,NP,NN) )        
        f.write("\n OR: %f, Fisher p_value: %f \t" %(oddsratio, pvalue))
        f.write("\n number PP: "+str(PP))
        f.write("\t\n PP:"+str(list(CHG_single.intersection(set(oncokb_list)))))
        
    gmt_dict["CHGs_in_oncoKB"]={}
    gmt_dict["CHGs_in_oncoKB"]["descriptor"]="baghera"
    gmt_dict["CHGs_in_oncoKB"]["genes"]=list(CHG_single.intersection(set(oncokb_list)))
        
    # all cancer list
    tot=len(all_genes_single)   
    PP=len(CHG_single.intersection(set(all_cancer_list)))
    PN=len(all_genes_single.intersection(set(all_cancer_list))  )-PP
    NP=len(CHG_single)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\t\n \n #### All cancer list enrichment")   
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(CHG_single.intersection(set(all_cancer_list)))))

    # oncoKB CGC
    tot=len(all_genes_single)   
    PP=len(CHG_single.intersection(set(CGC_list)))
    PN=len(all_genes_single.intersection(set(CGC_list))  ) -PP
    NP=len(CHG_single)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\t\n \n #### Sanger CGC enrichment")   
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(CHG_single.intersection(set(CGC_list)))))  
        
    # TSG enrichment
    totals=all_genes_single.intersection(set(oncokb_list))
    CHG_oncoKB=CHG_single.intersection(set(oncokb_list))
    
    
    tot=len(totals)   
    PP=len(CHG_oncoKB.intersection(set(TSG_list)))
    PN=len(all_genes_single.intersection(set(TSG_list))  ) -PP
    NP=len(CHG_oncoKB)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\n \n ### OncoKB TSG enrichment")   
        f.write("\t\n \n There are %d that are in oncokb and in our geneset" %tot)
        f.write("\t\n \n Of those %d are TSG" %(len(all_genes_single.intersection(set(TSG_list))  ) ))
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(CHG_oncoKB.intersection(set(TSG_list)))))  

    # oncogeni enrichment

    tot=len(totals)   
    PP=len(CHG_oncoKB.intersection(set(oncogene_list)))
    PN=len(all_genes_single.intersection(set(oncogene_list))  ) -PP
    NP=len(CHG_oncoKB)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\n \n ### OncoKB oncogene enrichment")   
        f.write("\t\n \n There are %d that are in oncokb and in our geneset" %tot)
        f.write("\t\n \n Of those %d are oncogene" %(len(all_genes_single.intersection(set(oncogene_list))  ) ))
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(CHG_oncoKB.intersection(set(oncogene_list)))))  
        

    ################# Enrichment for pancancer CHG ######################################
    panCHG_single=set([b for a in panCHG_table.name.str.replace(" ","").str.split(";") for b in a])
    
    gmt_dict["panCHGs"]={}
    gmt_dict["panCHGs"]["descriptor"]="baghera"
    gmt_dict["panCHGs"]["genes"]=list(panCHG_single)
    
    with open(stats_file, "a") as f:          
        f.write("\t\n \n ### Enrichment CHG")           
        f.write("\t\n \nThere are %d areas in total that are pancancer " %len(set(panCHG_table["name"].values)))  
        f.write("\t\n \nThere are %d genes in total that are pancancer " %len(panCHG_single) )      
        f.write("\t\n \n"+str(panCHG_single) )    
        


    # oncoKB annotated
    tot=len(all_genes_single)   
    PP=len(panCHG_single.intersection(set(oncokb_list)))
    PN=len(all_genes_single.intersection(set(oncokb_list))) -PP
    NP=len(panCHG_single)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False)     
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\t\n \n #### OncoKB annotated enrichment")   
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )        
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(panCHG_single.intersection(set(oncokb_list)))))
        
    # all cancer list
    tot=len(all_genes_single)   
    PP=len(panCHG_single.intersection(set(all_cancer_list)))
    PN=len(all_genes_single.intersection(set(all_cancer_list)))-PP
    NP=len(panCHG_single)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\t\n \n #### All cancer list enrichment----------------")   
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(panCHG_single.intersection(set(all_cancer_list)))))

    # oncoKB CGC
    tot=len(all_genes_single)   
    PP=len(panCHG_single.intersection(set(CGC_list)))
    PN=len(all_genes_single.intersection(set(CGC_list)))-PP
    NP=len(panCHG_single)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\t\n \n #### Sanger CGC enrichment----------------")   
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n PP:"+str(list(panCHG_single.intersection(set(CGC_list)))))  
        
    # TSG enrichment
    totals=all_genes_single.intersection(set(oncokb_list))
    panCHG_oncoKB=panCHG_single.intersection(set(oncokb_list))
    
    
    tot=len(totals)   
    PP=len(panCHG_oncoKB.intersection(set(TSG_list)))
    PN=len(all_genes_single.intersection(set(TSG_list))  ) -PP
    NP=len(panCHG_oncoKB)-PP
    NN=tot-PP-PN-NP
    
    #oddsratio, pvalue,_,_ = stats.chi2_contingency([[PP,NP],[PN,NN]],correction=False) 
    oddsratio, pvalue = stats.fisher_exact([[PP, PN], [NP, NN]], alternative="greater")               
    with open(stats_file, "a") as f:          
        f.write("\t\n \n #### OncoKB panCHGs TSG enrichment----------------")   
        f.write("\t\n \n There are %d that are in oncokb and in our geneset" %tot)
        f.write("\t\n \n Of those %d are TSG" %(len(all_genes_single.intersection(set(TSG_list))  ) ))
                
        
        f.write("\t\n Fisher, PP: %d, PN: %d, NP: %d, NN: %d" %(PP,PN,NP,NN) )
        f.write("\t\n OR: %f, Fisher p_value: %f" %(oddsratio, pvalue))
        
        f.write("\t\n number PP:"+str(PP))
        f.write("\t\n number oncogene and panCHG:"+str(len(panCHG_oncoKB.intersection(set(oncogene_list)))))
        f.write("\t\n"+str(panCHG_oncoKB.intersection(set(oncogene_list)) ))
        f.write("\t\n number oncokb annotated with unknown function and panCHG:"+str(len(panCHG_oncoKB.difference(set(oncogene_list)).difference(set(TSG_list)))))
        f.write("\t\n"+str(panCHG_oncoKB.difference(set(oncogene_list)).difference(set(TSG_list))) )
        
        f.write("\t\n PP:"+str(list(panCHG_oncoKB.intersection(set(TSG_list)))))  



def create_section1(group):

    """
    Code to create all table/plots for section 1
    """    

    print("plotting heritability")
    figure1.plot_Heritability(table_heritability_file, figs_output)
    print("plotting enrichment")   
    figure1.plot_Enrichment(group, figs_output)
    print("plotting minSNP")
    figure1.plot_minSNP(figs_output,tables_output+"SNPvsHG.csv")
    print("creating significant genes table")
    tab_significant=table_significant(group)        
    print("plotting prevalence")
    figure1.plot_HeritabilityPrevalence(group, tab_significant,table_prevalence_file, figs_output, tables_output)
    #print("plotting significant genes per cancer")    
    #figure1.plot_significant(figs_output, tab_significant, tables_output+'tab_liability.csv')

    
def create_section2(group):
    

    print("GO analysis")
    table_GO=GOOverlaps(group)
    figure2.plot_pathways_fisher(table_GO, figs_output,"GO",p_col="bh_pvalue",id_col="term",plotting_col="OR",threshold=0.1)
    figure2.pathways_fisher_volcano(table_GO, figs_output,"GO",p_col="bh_pvalue",id_col="term",plotting_col="OR",threshold=0.1)
    print("KEGG analysis")
    table_KEGG=KEGGOverlaps(group)
    figure2.plot_pathways_fisher(table_KEGG,figs_output,"KEGG",p_col="bh_pvalue",id_col="path",plotting_col="OR",threshold=0.1)
    print("oncoKB figures")
    #figure2.plot_oncoKB(tables_output+"oncoKB.csv",figs_output)
    
def create_section3(group):
    
    print("number of malignancies")
    panCancer_table=panCancerTable(group)
    figure3.plot_hist_mal_genes(figs_output, panCancer_table)

    pivot_binary=panCancer_binary(group)
    panCancer_bin=pivot_binary[1]
    panCancer_bin.to_csv(tables_output+'binary_panCancer.csv')

    col=len(panCancer_bin.columns)
    panCancer_bin["name"]=panCancer_bin.index

    panCancer_binary_annotated=annotateGeneTable(panCancer_bin)   
    panGenes_driver=panCancer_binary_annotated[panCancer_binary_annotated["OncoKB"]==True]

    panGenes_driver["OncoKB_sort"]=panGenes_driver["OncoKB"]+2*panGenes_driver["OncoKB_TSG"]+panGenes_driver["OncoKB_Oncogene"]
    panGenes_driver=panGenes_driver.sort_values(by=["OncoKB_sort"])

    panGenes_driver_cut=panGenes_driver.iloc[:,:col]
    print(panGenes_driver_cut)

    labels=[]
    mapping=["Unknown","Oncogene","TSG","TSG&OG"]
    for val in panGenes_driver["OncoKB_sort"]:
        labels.append(mapping[val-1])
    print(len(labels))
    print(len(panGenes_driver_cut))
    
    panGenes_driver_cut['Function']=labels



    panGenes_driver_cut.to_csv(tables_output+'binary_panCancer_oncoKB.csv')
    print("panCancer Enrichment")
    pivot=panCancerEnrichment(group)
    figure3.plot_h2percentage(figs_output,pivot[0],"CHGs")    
    figure3.plot_h2percentage(figs_output,pivot[1],"panCHGs")
    print("panCancer heatmap")
    panGenes=pivot[1]
    panGenes["name"]=panGenes.index
    panGenes_annotated=annotateGeneTable(panGenes)   
    panGenes_driver=panGenes_annotated[panGenes_annotated["OncoKB"]==True]
    panGenes_driver.to_csv(tables_output+'panCancer_oncoKB.csv')
    figure3.plot_driver_heatmap(figs_output, panGenes_driver)
    print("TSG Enrichment")
    TSG_list=cancerGenesList[cancerGenesList["OncoKB TSG"]=="Yes"]["Hugo Symbol"].values.tolist()
    pivot=createTablePivot(group, TSG_list, 'TSG')
    figure3.plot_h2percentage(figs_output,pivot[0],"TSG")    
    if len(pivot[1]>1):
        figure3.plot_h2percentage(figs_output,pivot[1],"TSG_CHGs")
    if len(pivot[2]>1):
        figure3.plot_h2percentage(figs_output,pivot[2],"TSG_panCHGs")
    print("Oncogene Enrichment")
    OG_list=cancerGenesList[cancerGenesList["OncoKB Oncogene"]=="Yes"]["Hugo Symbol"].values.tolist()
    pivot=createTablePivot(group, OG_list, 'OG')
    figure3.plot_h2percentage(figs_output,pivot[0],"OG")    
    if len(pivot[1]>1):
        figure3.plot_h2percentage(figs_output,pivot[1],"OG_CHGs")
    if len(pivot[2]>1):
        figure3.plot_h2percentage(figs_output,pivot[2],"OG_panCHGs")

        


#%%
###########################################################
###### GROUP ANALYSIS #####################################
###########################################################


def groupAnalysis():
    
    """ Function that analyses a group of malignancies calls all the functions to create tables and figures
    """
    
    #initialization of variables
    
    group=pd.DataFrame(columns=['Malignancy','name','StatsMax','P','bg_median','mi_median','Significant','weight'])   
    #correlations=pd.DataFrame(columns=['Malignancy','tau_corr','tau_pvalue','Wtau_corr','Wtau_pvalue','spearman_corr','spearman_pvalue','W_stats','W_pvalue'])   
    SNPvsHG=pd.DataFrame(columns=['Malignancy','sigSNP','HG','sigSNP_HG',"tot"])   
    oncoKB=pd.DataFrame(columns=['Malignancy','NC','oncogene_only','TSG_only','oncoKB_only',"oncogene_TSG","tot"])  
    annotated=pd.DataFrame()   

    names=set()
    
    for results in files:
        malType=results[len(folderInput)+len("genesResults_"):-4] #files suffix
        with open(results) as f:
            GENES = pd.read_table(f, sep=',')    
        
        GENES=GENES[["name",'SNPs','StatsMax',"P","bg_median","mi_median"]]
        GENES["weight"]=(GENES["bg_median"]-GENES["mi_median"])/GENES["mi_median"]
        GENES["Malignancy"]=malType
        GENES["Significant"]=GENES["P"]>0.99
        
        # all significant genes, splitted
        geneset=set([b for a in GENES[GENES["P"]>0.99].name.str.replace(" ","").str.split(";") for b in a])

        # creates gmt dictionary
        gmt_dict[malType]={}
        gmt_dict[malType]["descriptor"]="baghera"
        gmt_dict[malType]["genes"]=list(geneset)
        
        # appends the whole table to the group one
        group = group.append(GENES,ignore_index=True)
        
        # creates table with annotations from cancerGenesList
        annotatedGene=annotateGeneTable(GENES.copy())
        
        #annotatedGene.to_csv(folderOutput+"final_tables/genesResults_"+str(malType)+".csv",index=False)
        
        
        #SNPvsHG
        annotatedGene["Malignancy"]=malType
        annotatedGene["pvalue"]=1-stats.chi2.cdf(annotatedGene["StatsMax"], 1)
        sigSNP_list=annotatedGene[annotatedGene["pvalue"]<1e-8]["name"].values.tolist()
        HG_list=annotatedGene[annotatedGene["P"]>0.99]["name"].values.tolist()
        both=set(sigSNP_list).intersection(set(HG_list))
        
        SNPvsHG=SNPvsHG.append({'Malignancy': malType,
                            'sigSNP':len(sigSNP_list)-len(both),
                            'HG':len(HG_list)-len(both),
                            'sigSNP_HG':len(both),
                            "tot":len(sigSNP_list)+len(HG_list)-len(both),
                            'sigSNP_list':sigSNP_list , 
                            'NotBaghera': list(set(sigSNP_list).difference(set(HG_list))),
                            'sigSNP_HG_list':both}, ignore_index=True)

        #oncoKB, only significant genes, checks numbero of TSG/oncogene/both
        annotatedGene=annotatedGene[annotatedGene["Significant"]==True]          
        annotatedGene["both"]=annotatedGene["OncoKB_Oncogene"] & annotatedGene["OncoKB_TSG"] 
        
        oncogene_TSG=len( annotatedGene[annotatedGene["both"]==True] )
        oncogene_only=len(annotatedGene[(annotatedGene["OncoKB_Oncogene"]==True) & (annotatedGene["both"]==False)])
        TSG_only=len(annotatedGene[(annotatedGene["OncoKB_TSG"]==True) & (annotatedGene["both"]==False)]) 
        oncoKB_only=len( annotatedGene[(annotatedGene["OncoKB"]==True) & (annotatedGene["OncoKB_TSG"]==False) & (annotatedGene["OncoKB_Oncogene"]==False)] ) 
        NC=len(annotatedGene[annotatedGene["OncoKB"]==False] )
        
        oncoKB=oncoKB.append({'Malignancy': malType,'NC':NC,'oncogene_TSG':oncogene_TSG,'oncogene_only':oncogene_only,'TSG_only':TSG_only, "oncoKB_only":oncoKB_only,"tot":len(annotatedGene)}, ignore_index=True)
        
        annotated = annotated.append(annotatedGene,ignore_index=True)
        
        #GENES=GENES.sort_values(by=["P"])
        
        #a=GENES["StatsMax"].rank(ascending=False).tolist()
        #GENES["chi_rank"]=a.copy()
        #b=GENES["P"].rank(ascending=False).tolist()
        #GENES["P_rank"]=b.copy()        

        #tau_corr, tau_pvalue=stats.kendalltau(GENES["chi_rank"],GENES["P_rank"])
        #Wtau_corr, Wtau_pvalue=stats.weightedtau(GENES["chi_rank"],GENES["P_rank"])
        #GENESSig=GENES[GENES["P"]>0.99]
        #s_corr,s_pvalue=stats.spearmanr(GENESSig["chi_rank"],GENESSig["P_rank"], axis=0, nan_policy='propagate')
        #W_stats,W_pvalue=stats.wilcoxon(GENESSig["chi_rank"],GENESSig["P_rank"])
        #correlations = correlations.append({'Malignancy': malType,'tau_corr':tau_corr,'tau_pvalue':tau_pvalue,'Wtau_corr':Wtau_corr,'Wtau_pvalue':Wtau_pvalue,
        #                                    'spearman_corr':s_corr,'spearman_pvalue':s_pvalue,'W_stats':W_stats,'W_pvalue':W_pvalue}, ignore_index=True)
        
    #correlations.to_csv(folderOutput+"chi_P_correlations_99.csv", sep=",", mode='w')
    
    #plotBar_oncoKB(oncoKB.copy())
    #plotBarSNP_HG(SNPvsHG.copy()) 
    
    oncoKB.to_csv(tables_output+"oncoKB.csv",index=False)    
    SNPvsHG.to_csv(tables_output+"SNPvsHG.csv",index=False)
    annotated.to_csv(tables_output+"oncoKB_annotation.csv",index=False)
    figure2.plot_oncoKB(annotated,figs_output,tables_output)

    significant=group[group["P"]>0.99]

    significant.to_csv(tables_output+"CHG.csv",index=False)
    panCHGs=[name for name,g in group[group["P"]>0.99].groupby(["name"]) if len(g)>1]  
    panCHG_table=significant[significant["name"].isin(panCHGs)]
    panCHG_table.to_csv(tables_output+"pancancer_CHG.csv",index=False)


    with open(stats_file, "w") as f:          
        f.write("# Results statistics")   
    ######### VARIOUS STATS #############
    significant=group[group["P"]>0.99]
    bg_min= np.min(significant["bg_median"])
    bg_max= np.max(significant["bg_median"])
    bg_median= np.median(significant["bg_median"])
    
    weight_min= np.min(significant["weight"])
    weight_max= np.max(significant["weight"])
    weight_median= np.median(significant["weight"])
    
    with open(stats_file, "a") as f:          
        f.write("\n \n ### Cancer Heritability Genes statistics")       
        
        f.write("\n \n min gene heritability estimate: %f for %s in %s \t " %(bg_min,significant[significant["bg_median"]==bg_min]["name"].values[0],significant[significant["bg_median"]==bg_min]["Malignancy"].values[0]))
        
        f.write("\n max gene heritability estimate: %f %s in %s \t" %( bg_max,significant[significant["bg_median"]==bg_max]["name"].values[0],significant[significant["bg_median"]==bg_max]["Malignancy"].values[0]     ))
        f.write("\n median gene heritability estimate: %f \t" %(bg_median) )
        
        f.write("\n \n min fold enrichment: %f for %s in %s \t" %(weight_min,significant[significant["weight"]==weight_min]["name"].values[0],significant[significant["weight"]==weight_min]["Malignancy"].values[0] ))
        f.write("\n max fold enrichment: %f %s in %s \t" %( weight_max,significant[significant["weight"]==weight_max]["name"].values[0],significant[significant["weight"]==weight_max]["Malignancy"].values[0]    ))
        f.write("\n median foldenrichment: %f \t" %( weight_median) )
        
    
    
    oncoKB_enrichments(group, annotated, panCHG_table)
    DNArepairOverlaps(group, DNArepair_file)
    hallmarksOverlaps(group, cosmic_file)

    #if len(files)==2:
    #    figure3.compare_ranking(stats_file, tables_output, figs_output, group)


    create_section1(group)
    create_section2(group)
    create_section3(group)

def parserFunc():
    
    parser = argparse.ArgumentParser(description='Parse command line options.')
    parser.add_argument('--folderInput', '-i',  type=str,default='../../../results/ukbb_results/', help="Data Input, use the gene_results folder ")
    parser.add_argument('--filesOutput', '-o',  type=str,default='../../../results/final_results/', help="output folder for files results")
    parser.add_argument('--tablesOutput', '-to',  type=str,default='../../../results/final_results/tables/', help="output folder for table results")
    parser.add_argument('--figuresOutput', '-fo',  type=str,default='../../../results/final_results/figures/', help="output folder for figures results")
    parser.add_argument('--KEGGgmt', '-K', type=str, default='../../../data/external/PathwayAnalysis/c2.cp.kegg.v6.1.symbols.gmt' )
    parser.add_argument('--GOcsv', '-G', type=str, default='../../../data/external/PathwayAnalysis/goslim_human_16102019.csv' )
    parser.add_argument('--cancer_gene_list', '-C', type=str, default='../../../data/external/PathwayAnalysis/CancerGenesList.txt' )
    parser.add_argument('--actionable_variant', '-A', type=str, default='../../../data/external/PathwayAnalysis/allActionableVariants.txt' )
    parser.add_argument('--DNArepair', '-R', type=str, default='../../../data/external/PathwayAnalysis/DNArepair.csv' )
    parser.add_argument('--cosmic', '-c', type=str, default='../../../data/external/PathwayAnalysis/Census_Census_allTue Apr_10_11_11_42_2018.csv' )
    parser.add_argument('--table_heritability', '-H', type=str, default='../../../results/table_heritability.csv' )
    parser.add_argument('--table_prevalence', '-P', type=str, default='../../../data/external/prevalence.csv' )

    return parser


if __name__ == "__main__":
    
    """
    This script runs a complete post-hoc analysis on the results of BAGHERA.
    It creates all the tables and figures.
    """
    reset_matplotlib()

    parser = parserFunc()
    opts = parser.parse_args()
    
    #input files
    folderInput=opts.folderInput # folder with all the results files in csv (generated by baghera)
    files=g.glob(folderInput+"*.csv")
    table_heritability_file=opts.table_heritability
    table_prevalence_file=opts.table_prevalence
    DNArepair_file=opts.DNArepair
    cosmic_file=opts.cosmic

    #output folders
    files_output=opts.filesOutput #output folder for tables
    tables_output=opts.tablesOutput #output folder for tables
    figs_output=opts.figuresOutput #output folder for figures 



        
    with open(opts.actionable_variant,"r") as f:
        actionable=pd.read_table(f)
    with open(opts.cancer_gene_list,"rb") as f:       
        cancerGenesList= pd.read_table(f, sep='\t')        
    #TUMOR=TUMOR.drop_duplicates(subset=["Gene"],keep="first")

    #tab_prevalence_file='/home/viola/Documents/UniEssex/repos/project-cancer-heritability-gene/results/clean_results/group_results/prevalence.csv'
    #tab_prevalence_file='../results/clean_results/group_results/tab3_prevalence.csv'
   
    with open(opts.GOcsv,"r") as f:
        GO=pd.read_table(f,sep=",")
        
    KEGG=gp.parser.gsea_gmt_parser(opts.KEGGgmt)
    
    # Color settings
    color_blu=["#045a8d"]
    color_red=["#b2182b"]
    palette_binary = RdBu_4.mpl_colors[0::3]
    palette_sequential_blu = PuBu_8.mpl_colors #PuBu_9
    palette_sequential_red = Amp_20.mpl_colors #Reds_9
    #sns.palplot(palette_sequential_blu)

    stats_file=files_output+"stats.txt"
    
    gmt_dict={}
    groupAnalysis()
    printGMT(gmt_dict, files_output+"genes_lists.gmt")
    
    
