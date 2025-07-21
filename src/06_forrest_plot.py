import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import forestplot as fp
from datetime import datetime
import os



# Load the Excel file
def load_sheet_data(file_path):
    excel_data = pd.ExcelFile(file_path)
    sheet_data = {sheet: excel_data.parse(sheet) for sheet in excel_data.sheet_names}

    # Ensure column names are strings
    for sheet, df in sheet_data.items():
        df.columns = df.columns.map(str)

    return sheet_data, excel_data.sheet_names[0]


if __name__ == "__main__":
    

    today = datetime.now()   # Get date
    results_dir = "../results/"+today.strftime("%Y%m%d")  

    # Create the directory 'ihritik'
    try:
        os.makedirs(results_dir, exist_ok=True)
        print("Directory '%s' created successfully" % results_dir)
    except OSError as error:
        pass 

    cymap={
           "IL.1alpha":"IL-1$\\alpha$", 
           "IL.4":"IL-4", 
           "IL.5":"IL-5", 
           "IL.8":"IL-8", 
           "IL.10":"IL-10", 
           "IL.13":"IL-13", 
           "IL.17":"IL-17", 
           "IL.22":"IL-22", 
           "IL.24":"IL-24", 
           "IL.33":"IL-33", 
           "IL.1F7":"IL-37", 
           "IFN.g":"IFN-$\\gamma$", 
           "Eotaxin.3":"CCL-26", 
           "G.CSF":"GCSF", 
           "Periostin":"POSTN", 
           "SCGB1A.1":"SCGB1A1", 
           "TNF.alpha":"TNF-$\\alpha$", 
           "TSLP":"TSLP"
    }

    cellmap={'MS_DIFF_MONOS': "Monocytes", 'MS_DIFF_MAKROS':"Macrophages", 'MS_DIFF_NEUT':"Neutrophils",
             'MS_DIFF_EOS':"Eosinophils", 'MS_DIFF_LYM':"Lymphocytes", 'MS_DIFF_FLIMMEREPITHEL':"Ciliated epithelial cell"}

    
    groups = ["Severe", "Moderate", "Mild", "All"]
    grey_scale = ["#525252", "#969696", "#cccccc", "#d6d6d6"]
    sig_scale = ["#d7301f","#fc8d59","#fdcc8a", "#2b8cbe"]
    
    
        
    # Clinical Variable forrest plot
    
    fig = plt.figure()
    file_path = "../supplementary_tables/supplementary_table_{}.xlsx" 
    
    
    for i in range(3,12):
        print(i)
        sheet_data, name = load_sheet_data(file_path.format(i))
        print("processing {} file: {}".format(name, file_path.format(i)))
        data_dict={
        "All": sheet_data[name],
        "Mild": sheet_data['mild'],
        "Moderate": sheet_data['moderat'],
        "Severe": sheet_data['severe']
    }
    
        combined_data = []
        for group_name, df in data_dict.items():
                df_tmp = df.copy()
                df_tmp = df_tmp.iloc[range(16),:]
                df_tmp['Group'] = group_name
                combined_data.append(df_tmp)
        combined_data = pd.concat(combined_data, ignore_index=True)
        combined_data["Predictor"] = [cymap[c] for c in combined_data["Predictor"] ]
        cyts = combined_data["Predictor"].unique().tolist()[::-1]
        
        std_opts = dict(
            dataframe = combined_data,
            estimate = "Estimate",
            ll  ="CIlow", 
            hl = "CIhigh",
            varlabel="Predictor",
            model_col="Group",
        )
    
        
        _df, ax = fp.mforestplot(**std_opts,
                                 color_alt_rows=True,
                                 table=False,
                                 #right_annoteheaders=["Cytokine", "Severity Group"],
                                 xlabel="Coefficient (95% CI)",
                                 mcolor=["#2b8cbe","#fdcc8a","#fc8d59","#d7301f"][::-1],
                                 #xticks=[-1200,-600, 0, 600],
                                 return_df=True,
                                 despine=False,
                                 # Additional kwargs for customizations
                                 **{"marker": "D",  # set maker symbol as diamond
                                    "markersize": 35,  # adjust marker size
                                    "xlinestyle": (0, (10, 5)),  # long dash for x-reference line 
                                    "xlinecolor": "#808080",  # gray color for x-reference line
                                    "xtick_size": 12,  # adjust x-ticker fontsize
                                    "xlinestyle": (0, (10, 5)),  # long dash for x-reference line
                                    "xlinecolor": ".8",  # gray color for x-reference line
                                    "title":"{}".format(name),
                                   }                           
                                )
    
        # Access the PathCollection object from the Axes
        for i in range(4):
            
            scatter_collection = ax.collections[i]
            lines_collections = ax.collections[i+4]
    
            facecolors = scatter_collection.get_facecolor() 
    
            colors = []
            for c in cyts:
                #print(c, groups[i])
                #print(combined_data[(combined_data["Unnamed: 0"] == c) & (combined_data["Group"] == groups[i])]["Pr(>|z|)"])
                try:
                    p = float(combined_data[(combined_data["Predictor"] == c) & (combined_data["Group"] == groups[i])]["p_value)"])
                except:
                    p = float(combined_data[(combined_data["Predictor"] == c) & (combined_data["Group"] == groups[i])]["p_value"])
                #print("{}, {} pval={}".format(c,  groups[i], p))
                if p > 0.05:
                    colors.append(grey_scale[i])
                else:
                    colors.append(sig_scale[i])
                
    
            # Update the scatter plot with new sizes
            scatter_collection.set_facecolor(colors)
            lines_collections.set_color(colors)
    
        
        ax.set_title(name.replace("_", "/"))
        plt.tight_layout()
        plt.savefig(results_dir+'/figures/{}_forestplot.pdf'.format(name))

    # Compositional data Forrest plot
    file_path = "../supplementary_tables/supplementary_table_12.xlsx"
    df = pd.read_excel(file_path)
    df.cell_type.unique()

    # group by cell type 
    groups = df.groupby("cell_type")

    for cell_type, cdf in groups:
    
        df_tmp = cdf[cdf.predictor.isin(cymap.keys())]
        df_tmp["predictor"] = [cymap[c] for c in df_tmp["predictor"]]
    
        _df, ax = fp.forestplot(df_tmp,
                                 estimate= "beta_debiased",
                                 ll  ="ci_low", 
                                 hl = "ci_high",
                                 varlabel="predictor",
                                 color_alt_rows=True,
                                 table=False,
                                 #right_annoteheaders=["Cytokine", "Severity Group"],
                                 xlabel="Coefficient (95% CI)",
                                 #xticks=[-1200,-600, 0, 600],
                                 return_df=True,
                                 despine=False,
                                 ci_report=False,  # Turn off conf. int. reporting
                                 flush=False,  # Turn off left-flush of text
                                # Additional kwargs for customizations
                                 **{
                                    "figsize": [3,6],
                                    "marker": "x",  # set marker symbol as diamond
                                    "markersize": 36,  # adjust marker size
                                    "fontsize": 12,
                                    "xlinestyle": (0, (10, 5)),  # long dash for x-reference line 
                                    "xlinecolor": "#808080",  # gray color for x-reference line
                                    "xtick_size": 9.5,  # adjust x-ticker fontsize
                                    "xlinestyle": (0, (10, 5)),  # long dash for x-reference line
                                    "xlinecolor": ".8",  # gray color for x-reference line
                                    "title":"{}".format(cell_type),
                                    "col_spacing": 100
                                   }        
                                )
        # Access the PathCollection object from the Axes
        scatter_collection = ax.collections[1]
        lines_collections = ax.collections[0]
    
        facecolors = scatter_collection.get_facecolor() 
        colors =  [ grey_scale[-1] if r.p_value >= 0.05 else sig_scale[-1]
                    for j,r in df_tmp.iterrows() ]
    
        
        # Update the scatter plot with new sizes
        scatter_collection.set_facecolor(colors[::-1])
        lines_collections.set_color(colors[::-1])
        
        ax.set_title(cellmap[cell_type], fontdict={"fontsize":14})
        plt.tight_layout()
        plt.savefig(results_dir+'/figures/{}_forestplot.pdf'.format(cellmap[cell_type].replace(" ","_")))
