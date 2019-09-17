# -*- coding: utf-8 -*-

from jinja2 import Environment, FileSystemLoader
import os
import tornado.ioloop
import tornado.web
import pandas as pd
import numpy as np
#from scipy.optimize import fsolve
#from string import strip
from scripts.utility import make_responsive_plot, find_code, protein#, make_multi_plot
import json
import re
from bokeh.embed import components
from bokeh import layouts
import gc
#from bokeh.io import vform
#from bokeh.plotting import figure
#from bokeh.embed import components


#_WEB_SITE = 'http://localhost:8082/turnover'
_PORT = 8082
_PATH = os.path.dirname(os.path.abspath(__file__))
static_path=os.path.join(_PATH, 'static')

templateLoader = FileSystemLoader(searchpath="templates")
templateEnv = Environment(loader=templateLoader)

print 'load data'
#df = pd.DataFrame.from_csv('merge.csv')

bsf_data = pd.DataFrame.from_csv('Table_S2.csv')
bsf_data.columns = [n+'_leftdf' for n in  bsf_data.columns]
pcf_data = pd.DataFrame.from_csv('Table_S3.csv')
pcf_data.columns = [n+'_rightdf' for n in  pcf_data.columns]

def expand(in_df,col):
    for index, pg in zip(in_df.index.values,in_df[col]):
        pg_list = pg.split(';')
        if len(pg_list)>1:
            for prot in pg_list[1:]:
                if '__' in prot:
                    pass
                else:
                    temp = in_df.loc[index]
                    temp.name = prot
                    in_df = in_df.append(temp)
    return in_df

bsf_data = expand(bsf_data,'protein_groups_leftdf')
pcf_data = expand(pcf_data,'protein_groups_rightdf')

merge = bsf_data.join(pcf_data, how='outer')
merge['protein_group'] = merge['protein_groups_leftdf'].fillna('')+' '+merge['protein_groups_rightdf'].fillna('')
merge['protein_group_code'] = merge.protein_group.astype('category').cat.codes

merge.loc['Tb927.5.4480']

merge['index_bk']=merge.index.values
merge = merge.drop_duplicates(keep='first',subset='index_bk')
del merge['index_bk']

merge.loc['Tb927.5.4480']

id_dict = dict(zip(merge.index.values, merge['protein_group']))
print id_dict['Tb927.5.4480']


def grab(x):
    x=str(x).split(';')
    temp_id='none'
    temp_desc='none'
    for n in x:
        if 'ID=' in n:
            temp_id=n.replace('ID=','')
        if 'description=' in n:
            temp_desc=n.replace('description=','')
            temp_desc=temp_desc.replace('%2C',',')
    return temp_id, temp_desc
           
gff = pd.read_table('TriTrypDB-41_TbruceiTREU927.gff',header=None)
desc_dict = {'Tb427.BES40.22':'main vsg'}
for n in gff.iloc[:,-1]:
    temp = grab(n)
    desc_dict[temp[0]]=temp[1]
    

def get_table(merge):
    table = [] 
    for a,b,c,d in zip(merge['protein_group_code'], merge['new_half_life_leftdf'], merge['new_half_life_rightdf'], merge.index.values):
        #parse_text(df.loc[temp_id]['bs_exp_expected_double'])
        if d not in desc_dict:
            desc_dict[d]='none'

        item = dict(ProteinGroupId=d, 
                    Description = desc_dict[d],
                    BSF_r = b,
                    PCF_r = c,
                    ProteinGroupBk=a
                    )
        table.append(item)
    return table

    

table = get_table(merge)
print id_dict['Tb927.5.4480']
print table[0]
#id_dict = dict(zip(df['clean_ids'].values, df.index.values))
#bsf_avg =  pd.DataFrame.from_csv('BS/BS_deg_norm_avg.csv')
#bsf_std = pd.DataFrame.from_csv('BS/BS_deg_norm_std.csv')
#bsf_na = pd.DataFrame.from_csv('BS/BS_deg_norm_countna.csv')

#pcf_avg =  pd.DataFrame.from_csv('PC/PC_deg_norm_avg.csv')
#pcf_std = pd.DataFrame.from_csv('PC/PC_deg_norm_std.csv')
#pcf_na = pd.DataFrame.from_csv('PC/PC_deg_norm_countna.csv')

def func_def(x, N0, tau, offset):
    return N0 * np.exp( - x / tau ) + offset

#s1, s2 = make_multi_plot(df)
def make_out(request_prot):
    request_prot = request_prot.split(';')[0]
    temp_bsf_data = merge[ [n for n in merge.columns if n.endswith('_leftdf')] ]
    temp_bsf_data.columns = [n.replace('_leftdf','') for n in temp_bsf_data.columns]
    temp_bsf_data = temp_bsf_data.loc[request_prot]
    #print(temp_bsf_data)
    
    temp_pcf_data = merge[ [n for n in merge.columns if n.endswith('_rightdf')] ]
    temp_pcf_data.columns = [n.replace('_rightdf','') for n in temp_pcf_data.columns]
    temp_pcf_data = temp_pcf_data.loc[request_prot]
    #print(temp_pcf_data)
    #print(temp_pcf_data) 
    bsf_x_fit = np.array([0,0.5,1,2,4,8,12])#np.arange(0,13,1)
    popt = temp_bsf_data[['amplitude','tau','offset']]
    #print(request_prot, popt)
    if str(popt[0]) != 'nan':
        bsf_y_fit = func_def(bsf_x_fit, *popt)
        bsf_x_fit = ['x_fit']+[round(n,2) for n in list(bsf_x_fit)]
        bsf_y_fit = ['y_fit']+[round(n,2) if n < 1 else 1 for n in list(bsf_y_fit)]
    else:
        bsf_x_fit = ['x_fit','']
        bsf_y_fit = ['y_fit','']
    
    pcf_x_fit = np.array([0, 0.25, 0.5, 1, 2, 4, 8, 20, 28])#np.arange(0,29,1)
    popt = temp_pcf_data[['amplitude', 'tau', 'offset']]
    if str(popt[0]) != 'nan':
        pcf_y_fit = func_def(pcf_x_fit, *popt)
        pcf_x_fit = ['x_fit']+[round(n,2) for n in list(pcf_x_fit)]
        
        pcf_y_fit = ['y_fit']+[round(n,2) if n < 1 else 1 for n in list(pcf_y_fit)] 
    else:
        pcf_x_fit=['x_fit','']
        pcf_y_fit=['y_fit','']
        
    html_df = pd.DataFrame()
    html_df['Fit Values'] = ['Amplitude', 'Tau', 'Offset', 'Half Life', 'RMSE' ]
    #print(temp_bsf_data[['N0','tau','offset','Exponential_Half_Life']])
    html_df['BSF'] = temp_bsf_data[['amplitude','tau','offset', 'new_half_life', 'exp_rmse']].values
    html_df['PCF'] = temp_pcf_data[['amplitude','tau','offset', 'new_half_life', 'exp_rmse']].values
    
    html_df = html_df.round(2)
    html_table = html_df.to_html(index=False, classes='display')
    html_table = re.sub(' display', ' display" id="fit_table"', html_table)
    
    start_data = 13
    #print temp_bsf_data.values[0+start_data:7+start_data]
    out = {
            
            'bsf':{
            'request_prot':request_prot,
            'data_1':['data_1']+[ round(n,2) if str(n) != 'nan' else '' for n in temp_bsf_data.values[0+start_data:7+start_data]],
            'data_2':['data_2']+[ round(n,2) if str(n) != 'nan' else '' for n in temp_bsf_data.values[7+start_data:14+start_data]],
            'data_3':['data_3']+[ round(n,2) if str(n) != 'nan' else '' for n in temp_bsf_data.values[14+start_data:21+start_data]],                
            'x_point':['x_point']+[0,0.5,1,2,4,8,12],
            'x_fit':bsf_x_fit,                        
            'y_fit':bsf_y_fit,
                    },
                    
            'pcf':{
            'request_prot':request_prot,
            'data_1':['data_1']+[ round(n,2) if str(n) != 'nan' else '' for n in temp_pcf_data.values[0+start_data:9+start_data] ],
            'data_2':['data_2']+[ round(n,2) if str(n) != 'nan' else '' for n in temp_pcf_data.values[9+start_data:18+start_data] ],
            'data_3':['data_3']+[ round(n,2) if str(n) != 'nan' else '' for n in temp_pcf_data.values[18+start_data:27+start_data] ],                
            'x_point':['x_point']+[0,0.25,0.5,1,2,4,8,20,28],
            'x_fit':pcf_x_fit,                        
            'y_fit':pcf_y_fit,
                    }, 
                    
            'desc':desc_dict[request_prot],
            'fit_table':html_table   
           }
    return out
    
print 'load big data'    
all_data = {}

for item in table:
    #print item
    temp_id = item['ProteinGroupId']
    out = make_out(temp_id)
    del out['fit_table']
    all_data[temp_id]=out
    #break
    
print 'start server'  
print all_data['Tb05.5K5.110']
print all_data['Tb927.5.4460']  
#print df.loc['Tb05.5K5.210;Tb927.5.4570']
       
class MainHandler(tornado.web.RequestHandler):
    def get(self):
        template_file = "turnover.html"
        template = templateEnv.get_template(template_file)
        output = template.render(
                title="welcome to protein turnover", 
                table=table
                )
        self.write(output)


class GetData(tornado.web.RequestHandler):
    def get(self):
        request_prot = self.request.uri.split('/')[-1]
        out = make_out(request_prot)
        self.write(json.dumps(out))
        
class ProteinPlot(tornado.web.RequestHandler):
    
    def get(self):
        start_data=13
        request_prot = self.request.uri.split('/')[-1] 
        request_prot = request_prot.split(';')[0]
        
        
        temp_bsf_data = merge[ [n for n in merge.columns if n.endswith('_leftdf')] ]
        temp_bsf_data.columns = [n.replace('_leftdf','') for n in temp_bsf_data.columns]
        temp_bsf_data = temp_bsf_data.loc[request_prot]
        
        temp_pcf_data = merge[ [n for n in merge.columns if n.endswith('_rightdf')] ]
        temp_pcf_data.columns = [n.replace('_rightdf','') for n in temp_pcf_data.columns]
        temp_pcf_data = temp_pcf_data.loc[request_prot]

        prot_id = request_prot
        desc = desc_dict[request_prot]
        request_tag = self.request.uri.split('/')[-2] 
        
        print prot_id, desc, request_prot
        
        if request_tag == 'BSF':
            Ycode=find_code(temp_bsf_data.isnull()[0+start_data:21+start_data])
            N0=temp_bsf_data['amplitude']
            tau=temp_bsf_data['tau']
            offset=temp_bsf_data['offset'] 
            half_life=temp_bsf_data['new_half_life'] 
            r_square=0#temp_bsf_data['exp_r2'] 
            X=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]

            bsf = pd.DataFrame()
            bsf['A']=temp_bsf_data.values[0+start_data:7+start_data]
            bsf['B']=temp_bsf_data.values[7+start_data:14+start_data]
            bsf['C']=temp_bsf_data.values[14+start_data:21+start_data]
            bsf['mean']=bsf.iloc[:,0:3].mean(axis=1,skipna=True)
            bsf['std']=bsf.iloc[:,0:3].std(axis=1,skipna=True)

            Y=bsf['mean'].values
            Yerr=bsf['std'].values
            prot = protein(prot_id, desc , N0 , tau, offset, half_life, r_square, X, Y, Yerr, Ycode)
            prot.replace_na()
            
        if request_tag == 'PCF':
            Ycode=find_code(temp_pcf_data.isnull().values[0+start_data:27+start_data])
            N0=temp_pcf_data['amplitude']
            tau=temp_pcf_data['tau']
            offset=temp_pcf_data['offset'] 
            half_life=temp_pcf_data['new_half_life'] 
            r_square=0#temp_pcf_data['exp_r2'] 
            X=[0,0.25,0.5,1,2,4,8,20,28]
            
            pcf = pd.DataFrame()
            pcf['A']=temp_pcf_data.values[0+start_data:9+start_data]
            pcf['B']=temp_pcf_data.values[9+start_data:18+start_data]
            pcf['C']=temp_pcf_data.values[18+start_data:27+start_data]
            pcf['mean']=pcf.iloc[:,0:3].mean(axis=1,skipna=True)
            pcf['std']=pcf.iloc[:,0:3].std(axis=1,skipna=True)


            Y=pcf['mean'].values
            Yerr=pcf['std'].values
            Ycode=find_code(temp_pcf_data.isnull().values[0+start_data:27+start_data])
            prot = protein(prot_id, desc, N0, tau, offset,half_life,r_square,X,Y,Yerr,Ycode)
            prot.replace_na()            
            
        
        #make_responsive_plot(prot, tag='BSF')        
        
        script, div = make_responsive_plot(prot, tag=request_tag)
        
        template_file = "protein_plot.html"
        template = templateEnv.get_template(template_file)
        output = template.render(
                script=script, 
                div=div,
                protein_id=request_prot,
                protein_desc=desc,
                start_hl=tau/np.log(2),
                start_r=half_life
                )
        self.write(output)        
        

class MultiPlot(tornado.web.RequestHandler):
    def get(self):
        #
        template_file = "multi_plot.html"
        template = templateEnv.get_template(template_file)
        
        output = template.render(
                title="welcome to protein turnover", 
                table=table,
                all_protein_data = json.dumps(all_data),
                )
        
        self.write(output)        
        #base_script, base_div = components(layouts.row(s1, s2))  
        
        
        #output = template.render(script=base_script, div=base_div)
        #output = 'work in progess'
        #self.write(output)   
        
        
handlers=[
    (r"/static/(.*)", tornado.web.StaticFileHandler, {'path':static_path} ),
    (r"/turnover", MainHandler),
    (r"/get_data/.*", GetData),
    (r"/protein_plot/.*", ProteinPlot),
    (r"/multi_plot", MultiPlot),    
    ]


settings = {'debug': True, }


def make_app():
    return tornado.web.Application(handlers , **settings) 
    

if __name__ == "__main__":
    #pass
    # Setup the server
    app = make_app()
    app.listen(_PORT)
    tornado.ioloop.IOLoop.instance().start()
    #pass