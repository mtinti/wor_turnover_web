from bokeh import layouts #import layout
from bokeh.io import show
from bokeh.plotting import figure
import bokeh.models as bkm
from bokeh.embed import components
import numpy as np
#from bokeh.models import Range1d
import pandas as pd
from itertools import izip
from bokeh.models import Legend
from bokeh.models.glyphs import Line
from bokeh.models.glyphs import MultiLine

def model_deg(x, N0, tau, c):
    return N0 * np.exp( -x/tau ) + c

def model_syn(x, N0, tau, c):
    return 1 - ( N0 * np.exp( -x/tau ) +c )

def find_code(X):
    temp = pd.DataFrame()
    limit  = int(X.shape[0]/3.0)
    print list(X)
    print limit
    temp['A']=X[0:limit]
    temp['B']=X[limit:limit*2]
    temp['C']=X[limit*2:limit*3]

    res = temp.sum(axis=1) 
    print 'res',res
    return 3-res
    
class protein():
    def __init__(self, prot_id,
                 desc, N0, 
                 tau, 
                 offset, 
                 half_life, 
                 r_square, 
                 X, Y, Yerr, Ycode):
        
        self.id = prot_id
        self.desc  = desc
        self.N0 = N0
        self.tau = tau
        self.offset = offset
        self.half_life = half_life
        self.r_square = r_square
        self.X = X
        self.Y = Y
        self.Yerr = Yerr
        self.Ycode = Ycode
    
    def replace_na(self):
        if str(self.N0) == 'nan':
            self.N0=0
        if str(self.tau) == 'nan':
            self.tau=0 
        if str(self.offset) == 'nan':
            self.offset=0  
        if str(self.r_square) == 'nan':
            self.r_square=0  
        if str(self.half_life) == 'nan':
            self.half_life=0  
        self.X  = np.array([n  if str(n) != 'nan' else 0 for n in self.X ] )
        self.Y  = np.array([n  if str(n) != 'nan' else 0 for n in self.Y ] )
        self.Yerr  = np.array([n  if str(n) != 'nan' else 0 for n in self.Yerr ] )
        self.Ycode  = np.array([n  if str(n) != 'nan' else 0 for n in self.Ycode ]  )         
        #pass
    def __str__(self):
        
        data = [self.id, 
               self.desc, 
               self.N0 ,
               self.tau ,
               self.offset,
               self.half_life, 
               self.r_square, 
               self.X,
               self.Y,
               self.Yerr,
               self.Ycode
        ]
        label = ['id', 
               'desc', 
               'N0',
               'tau',
               'offset',
               'half_life', 
               'r_square', 
               'X',
               'Y',
               'Yerr',
               'Ycode'
        ]
        
        res = [str(a)+': '+str(b) for a,b in zip(label,data)]
        return '\n'.join(res)

def get_base_plot(prot, title='syn and deg curves'):
    #print temp_prot.name
    #print temp_prot.description
    fig_fit = figure(width=600, height=400,  title=title, 
                      toolbar_location="above")
    circle_data_deg_x = []
    circle_data_deg_y = []
    diamond_data_deg_x = []
    diamond_data_deg_y = []
    square_data_deg_x = []
    square_data_deg_y = []    
    for x,y,code in zip(prot.X,
                        prot.Y,
                        prot.Ycode):
        if code == 3:
            circle_data_deg_x.append(x)
            circle_data_deg_y.append(y)
        elif code == 2:
            diamond_data_deg_x.append(x)
            diamond_data_deg_y.append(y)                  
        elif code == 1:
            square_data_deg_x.append(x)
            square_data_deg_y.append(y)
    
    source_circle_deg = bkm.ColumnDataSource(data=dict(x=circle_data_deg_x, 
                                                y=circle_data_deg_y))
    source_diamond_deg = bkm.ColumnDataSource(data=dict(x=diamond_data_deg_x, 
                                                y=diamond_data_deg_y))
    source_square_deg = bkm.ColumnDataSource(data=dict(x=square_data_deg_x, 
                                                y=square_data_deg_y))                     
    
    fig_fit.scatter(x='x',y='y', source=source_square_deg, marker="square", color='red', size=0, fill_alpha=0.2)
    square = bkm.Square(x='x', y='y', fill_color='red', size=8, fill_alpha=0.2)
    square_r = fig_fit.add_glyph(source_or_glyph=source_square_deg, glyph=square) 
    square_hover = bkm.HoverTool(renderers=[square_r], tooltips=[('time', '@x'), ('N', '@y')])    
    fig_fit.add_tools(square_hover)

    fig_fit.scatter(x='x',y='y', source=source_diamond_deg, marker="diamond", color='red', size=0, fill_alpha=0.2)
    diamond = bkm.Diamond(x='x', y='y', fill_color='red', size=12, fill_alpha=0.2 )
    diamond_r = fig_fit.add_glyph(source_or_glyph=source_diamond_deg, glyph=diamond)                                    
    diamond_hover = bkm.HoverTool(renderers=[diamond_r],tooltips=[('time', '@x'), ('N', '@y')])
    fig_fit.add_tools(diamond_hover)
            
    fig_fit.scatter(x='x',y='y', source=source_circle_deg, marker="circle", color='red', size=0, fill_alpha=0.2)
    circle = bkm.Circle(x='x', y='y', fill_color='red', size=8, fill_alpha=0.2 )
    circle_r = fig_fit.add_glyph(source_or_glyph=source_circle_deg, glyph=circle)                                    
    circle_hover = bkm.HoverTool(renderers=[circle_r], tooltips=[('time', '@x'), ('N', '@y')])
    fig_fit.add_tools(circle_hover)
    
    
    circle_data_syn_x = []
    circle_data_syn_y = []
    diamond_data_syn_x = []
    diamond_data_syn_y = []
    square_data_syn_x = []
    square_data_syn_y = []            
    complement_Y = 1-prot.Y
    for x,y,code in zip(prot.X,
                        complement_Y,
                        prot.Ycode):
                            
        if code == 3:
            circle_data_syn_x.append(x)
            circle_data_syn_y.append(y)                    
            #s2.scatter(x,y,marker="circle",color='green',size=8)
        elif code == 2:
            diamond_data_syn_x.append(x)
            diamond_data_syn_y.append(y)                     
            #s2.scatter(x,y,marker="diamond",color='green',size=12, fill_alpha=0.2)
        elif code == 1:
            square_data_syn_x.append(x)
            square_data_syn_y.append(y)                    
            #s2.scatter(x,y,marker="square",color='green',size=8, fill_alpha=0.2 )        
    
    source_circle_syn = bkm.ColumnDataSource(data=dict(x=circle_data_syn_x, 
                                        y=circle_data_syn_y))
    source_diamond_syn = bkm.ColumnDataSource(data=dict(x=diamond_data_syn_x, 
                                        y=diamond_data_syn_y))
    source_square_syn = bkm.ColumnDataSource(data=dict(x=square_data_syn_x, 
                                        y=square_data_syn_y))

    fig_fit.scatter(x='x',y='y', source=source_square_syn, marker="square", color='green', size=0, fill_alpha=0.2)
    square = bkm.Square(x='x', y='y', fill_color='green', size=8, fill_alpha=0.2)
    square_r = fig_fit.add_glyph(source_or_glyph=source_square_syn, glyph=square)                                    
    square_hover = bkm.HoverTool(renderers=[square_r], tooltips=[('time', '@x'), ('N', '@y')])
    fig_fit.add_tools(square_hover)

    fig_fit.scatter(x='x',y='y', source=source_diamond_syn, marker="diamond", color='green', size=0, fill_alpha=0.2)
    diamond = bkm.Diamond(x='x', y='y', fill_color='green', size=12, fill_alpha=0.2 )
    diamond_r = fig_fit.add_glyph(source_or_glyph=source_diamond_syn, glyph=diamond)                                    
    diamond_hover = bkm.HoverTool(renderers=[diamond_r], tooltips=[('time', '@x'), ('N', '@y')])
    fig_fit.add_tools(diamond_hover) 

    fig_fit.scatter(x='x',y='y', source=source_circle_syn, marker="circle", color='green', size=0, fill_alpha=0.2)           
    circle = bkm.Circle(x='x', y='y', fill_color='green', size=8, fill_alpha=0.2 )       
    circle_r = fig_fit.add_glyph(source_or_glyph=source_circle_syn, glyph=circle)                                    
    circle_hover = bkm.HoverTool(renderers=[circle_r], tooltips=[('time', '@x'), ('N', '@y')])
    fig_fit.add_tools(circle_hover)            
    
    err_y = []
    err_x = []
    for x, y, yerr in zip(prot.X, prot.Y, prot.Yerr):
        err_x.append((x, x))
        err_y.append((y - yerr, y + yerr))
    fig_fit.multi_line(err_x , err_y,  color='red', line_dash="4 4")
    
    err_y = []
    err_x = []
    for x, y, yerr in zip(prot.X, complement_Y, prot.Yerr):
        err_x.append((x, x))
        err_y.append((y - yerr, y + yerr))    
    fig_fit.multi_line(err_x, err_y,  color='green',  line_dash="4 4") 
    return fig_fit



def make_responsive_plot(prot, tag):
    #output_file("test_deg_def.html", title="line deg")
    title = tag
        
    fig_fit = get_base_plot(prot, title)
    #fig_fit.set(y_range=bkm.Range1d(-0.1, 1.1))
    x_lim = prot.X[-1]+0.1
    X = np.arange(-0.1,x_lim,0.1)
    
    
    
    Ys = model_syn(X,prot.N0,prot.tau,prot.offset)
    Yd = model_deg(X,prot.N0,prot.tau,prot.offset)
    
    params_start = [prot.N0, prot.tau, prot.offset, title]
    params = [prot.N0, prot.tau, prot.offset, title] 

    diff = len(Yd)-len(params)
    add = [np.nan]*diff
    params_start+=add
    params+=add
    
    
    
    
    source = bkm.ColumnDataSource(data=dict(X = X, Ys = Ys,
                                        Yd = Yd,
                                        params_start = params_start,
                                        params = params,
                                        ) )
    
    code="""
                var data = source.get('data');
                var params_start = data['params_start']
                var X = data['X'];
                var Ys = data['Ys'];
                var Yd = data['Yd'];
                var params = data['params'];

                var title = params[3];
   
                var f = cb_obj.get('value');

                if (cb_obj.get('title') == 'N0') {N0 = cb_obj.get('value'); params[0] = N0;};
                if (cb_obj.get('title') == 'tau') {tau = cb_obj.get('value'); params[1] = tau; };
                if (cb_obj.get('title') == 'offset') {offset = cb_obj.get('value'); params[2] = offset;};
                
                var N0 = parseFloat(params[0]);
                var tau = parseFloat(params[1]);
                var offset = parseFloat(params[2]);
                
                //console.log(params, 'here parms');
                //console.log(N0, 'here N0');
                //console.log(tau, 'here tau');
                //console.log(offset, 'here offset');
                
                
                //console.log(Ys, 'here Ys before');
                //console.log(Yd, 'here Yd before');                
                
                for (i = 0; i < X.length; i++) {
                    Ys[i] =  1 - ((N0 * Math.exp( -X[i]/tau )) +  offset) ;
                }
                 
                for (i = 0; i < X.length; i++) {
                    Yd[i] =  (N0 * Math.exp( -X[i]/tau )) + offset  ;
                }        
                
                
                var last_element = X[X.length - 1];
                var residual =  (N0 * Math.exp( -last_element/tau )) + offset 
                
                
                //console.log(Ys, 'here Ys after');
                //console.log(Yd, 'here Yd after');
                
                
                source.trigger('change');
                console.log(N0, 'here N0 after');
                console.log(tau, 'here tau after'); 
                console.log(offset, 'here offset after'); 
                console.log( tau / Math.log(2) ,'here half life')
                //$("field_name_1").update(tau/Math.log(2));
                //$("field_name_2").update(tau);
                //$("field_name_3").update(offset);
                
                
                $('#field_name_1').text(Math.round(tau/Math.log(2)* 100)/100);
                $('#field_name_2').text(Math.round(residual * 100)/100);
                //$('#field_name_3').text(offset);
                
                //document.getElementById("field_name_1").innerHTML = N0;
                
                 // $( ".hello" ).remove();
                //var $new = $('<div class="hello"></div>').appendTo('.insert_body');
                //http://stackoverflow.com/questions/247483/http-get-request-in-javascript
                
                //function httpGet(theUrl)
                //{
                    //var xmlHttp = new XMLHttpRequest();
                    //xmlHttp.open( "GET", theUrl, false ); // false for synchronous request
                    //xmlHttp.send( null );
                    //return xmlHttp.responseText;
                //}
                //var address = '{web_address}find_intersection?N0='+N0+'&tau='+tau+'&c='+offset
                //var res = httpGet(address);
                //console.log(res, 'test');
                
                //if (title == 'Bloodstrem'){
                        //console.log('new_intersection_bs');
                        //$("#new_intersection_bs").html('new intersection: '+res);                
                //}
                //if (title == 'Procyclic'){
                       // $("#new_intersection_pc").html('new intersection: '+res);                
                //}
                
                
                //res = httpGet(address)
                //console.log(res, 'test');
                
                //ras = $.get("address");
                //req.open('GET', address, false);    
                //
                //console.log(address);
                """


    #code = code.replace('{web_address}', web_adress)
    callback = bkm.CustomJS(args=dict(source=source), code=code)     
    
    
    sliderN0 = bkm.Slider(start=0, end=1, value=prot.N0, step=.01, title="N0",  name = 'N0', callback=callback)
    sliderTau = bkm.Slider(start=0.1, end=50, value=prot.tau, step=.01, title="tau",  name = 'tau', callback=callback)
    sliderOffSet = bkm.Slider(start=0, end=1, value=prot.offset, step=.01, title="offset",  name = 'offset', callback=callback)
    syn = fig_fit.line(x='X', y='Ys', source=source, color='green')
    deg = fig_fit.line(x='X', y='Yd',  source=source, color='red')  

    legend = Legend(items=[ ("syn", [syn]), ("deg" , [deg])], location=(10, -30))
    #fig_fit.legend.location = (0, -30)
    #layout = vform(fig_fit)
    #script, div = components(layout)
    #layout =  vform()
    #script_bar, div_bar = components(layout)
    
    #show(layout(
            #[[fig_fit],[sliderN0],[sliderTau],[sliderOffSet]]

            #))
    fig_fit.add_layout(legend, 'right')
    #show(layouts.column(fig_fit, sliderN0, sliderTau, sliderOffSet))
    script, div = components(layouts.column(fig_fit, sliderN0, sliderTau, sliderOffSet))
    return script, div
    #print 'done'
    #return script, div, script_bar, div_bar





def make_multi_plot(merge):
    #s1 = figure(plot_width=400, plot_height=400, title=None)
    #s2 = figure(plot_width=400, plot_height=400, title=None)

    #merge = pd.DataFrame.from_csv('merge.csv')
    condition_1 = ((merge['bs_exp_offset']>=0) & (merge['bs_exp_offset']<=1))
    condition_2 = ((merge['bs_exp_rsquared']>=0.8) & (merge['bs_exp_rsquared']<=1))
    condition_3 = ((merge['select_line_aic_bs'] == False) & (merge['select_line_r_bs']== False))
    bs = merge[condition_1 & condition_2 & condition_3]
    X_bs = np.arange(-0.1, 13, 0.5) 
    source = bkm.ColumnDataSource(
            dict(
        xs=[X_bs for i in bs.index.values],
        ys=[model_deg(X_bs, bs.loc[n]['bs_exp_N0'], bs.loc[n]['bs_exp_tau'], bs.loc[n]['bs_exp_offset']) for n in bs.index.values],
            )
        )
    
    xdr = bkm.DataRange1d()
    ydr = bkm.DataRange1d()    

    s1 = bkm.Plot(title=None, x_range=xdr, y_range=ydr, plot_width=500, plot_height=500, toolbar_location='below')
    glyph = MultiLine(xs="xs", ys="ys", line_color="gray", line_width=1, line_alpha=0.1)
    s1.add_glyph(source, glyph)  

    yaxis = bkm.LinearAxis()
    xaxis = bkm.LinearAxis()
    s1.add_layout(xaxis, 'below')
    s1.add_layout(yaxis, 'left')                
                      
    condition_1 = ((merge['pc_exp_offset']>=0) & (merge['pc_exp_offset']<=1))
    condition_2 = ((merge['pc_exp_rsquared']>=0.8) & (merge['pc_exp_rsquared']<=1))
    condition_3 = ((merge['select_line_aic_pc'] == False) & (merge['select_line_r_pc']==False))
    pc = merge[condition_1 & condition_2 & condition_3]
    X_pc = np.arange(-0.1, 28, 0.5)
    source2 = bkm.ColumnDataSource(
            dict(
        xs=[X_pc for i in pc.index.values],
        ys=[model_deg(X_pc, pc.loc[n]['pc_exp_N0'], pc.loc[n]['pc_exp_tau'], pc.loc[n]['pc_exp_offset']) for n in pc.index.values],
            )
        )
    
    xdr = bkm.DataRange1d()
    ydr = bkm.DataRange1d()        
    s2 = bkm.Plot(title=None, x_range=xdr, y_range=ydr, plot_width=500, plot_height=500, toolbar_location='below')
    glyph = MultiLine(xs="xs", ys="ys", line_color="gray", line_width=1, line_alpha=0.1)
    
    s2.add_glyph(source2, glyph)  
    
    yaxis = bkm.LinearAxis()
    xaxis = bkm.LinearAxis()
    s2.add_layout(xaxis, 'below')
    s2.add_layout(yaxis, 'left') 
    
    #s2.add_layout(xaxis, 'below')
    #s2.add_layout(yaxis, 'left')
    
    
    
    
    s1.add_layout(bkm.Grid(dimension=0, ticker=xaxis.ticker))
    s1.add_layout(bkm.Grid(dimension=1, ticker=yaxis.ticker))    
    s2.add_layout(bkm.Grid(dimension=0, ticker=xaxis.ticker))
    s2.add_layout(bkm.Grid(dimension=1, ticker=yaxis.ticker))
    
    
    return s1, s2

if __name__ == '__main__':
    #X = np.arange(-0.1, 28, 0.5)
    #s2.multi_line([X for n in pc.index.values],
                  #[model_deg(X, pc.loc[n]['pc_exp_N0'], pc.loc[n]['pc_exp_tau'], pc.loc[n]['pc_exp_offset']) for n in pc.index.values],
             #color="gray", alpha=0.1, line_width=1)
    #script, div = components(layouts.row(s1, s2))
    #s1.multi_line([X for n in bs.index.values],
                  #[model_deg(X, bs.loc[n]['bs_exp_N0'], bs.loc[n]['bs_exp_tau'], bs.loc[n]['bs_exp_offset']) for n in bs.index.values],
             #color="gray", alpha=0.1, line_width=1)  
    #merge = pd.DataFrame.from_csv('../merge.csv')
    #make_multi_plot()
     
    def test_prot(request_prot = 'Tb927.5.4520'): 
        bsf_data = pd.DataFrame.from_csv('../Table_S2.csv')
        pcf_data = pd.DataFrame.from_csv('../Table_S3.csv')

        id_dict = {}
        desc_dict = {}
        for a,b,c in zip(bsf_data.index.values, bsf_data['protein_groups'], bsf_data['desc']):
            for prot in b.split(';'):
                id_dict[prot]=a
                desc_dict[prot]=c
        for a,b,c in zip(pcf_data.index.values, pcf_data['protein_groups'], pcf_data['desc']):
            for prot in b.split(';'):
                if prot not in id_dict:
                    id_dict[prot]=a
                    desc_dict[prot]=c  

        
        
        desc = desc_dict[request_prot]

        if id_dict[request_prot] in bsf_data.index.values:
            temp_bsf_data = bsf_data.loc[id_dict[request_prot]]
        else:
            temp_bsf_data = bsf_data.loc[bsf_data.index.values[0]]
            temp_bsf_data._update_inplace=np.nan
    
        if id_dict[request_prot] in pcf_data.index.values:
            temp_pcf_data = pcf_data.loc[id_dict[request_prot]]
        else:
            temp_pcf_data = pcf_data.loc[bsf_data.index.values[0]]
            temp_pcf_data._update_inplace=np.nan

          
        #merge = pd.DataFrame.from_csv('../merge.csv')
        #bsf_avg =  pd.DataFrame.from_csv('../BS/BS_deg_norm_avg.csv')
        #bsf_std = pd.DataFrame.from_csv('../BS/BS_deg_norm_std.csv')
        #bsf_na = pd.DataFrame.from_csv('../BS/BS_deg_norm_countna.csv')
        #prot_id ='Tb927.10.1070:mRNA-p1'

        
        Ycode=find_code(temp_bsf_data.isnull().values[0:21])
        N0=temp_bsf_data['amplitude']
        tau=temp_bsf_data['tau']
        offset=temp_bsf_data['offset'] 
        half_life=temp_bsf_data['exp_residual_life'] 
        r_square=temp_bsf_data['exp_r2'] 
        X=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]

        bsf = pd.DataFrame()
        bsf['A']=temp_bsf_data.values[0:7]
        bsf['B']=temp_bsf_data.values[7:14]
        bsf['C']=temp_bsf_data.values[14:21]
        bsf['mean']=bsf.iloc[:,0:3].mean(axis=1,skipna=True)
        bsf['std']=bsf.iloc[:,0:3].std(axis=1,skipna=True)

        Y=bsf['mean'].values
        Yerr=bsf['std'].values
        prot = protein(request_prot, desc , N0 , tau, offset, half_life, r_square, X, Y, Yerr, Ycode)
        prot.replace_na()
        print(prot)
    
    
    
        script, div = make_responsive_plot(prot, tag='BSF')
        return script, div
        #script, div = make_responsive_plot(prot, tag=request_tag)
    
    request_prot =    'Tb927.5.4520' 
    script, div = test_prot(request_prot = request_prot)
    
    from jinja2 import Environment, FileSystemLoader
    templateLoader = FileSystemLoader(searchpath="../templates")
    templateEnv = Environment(loader=templateLoader)

    template_file = "protein_plot.html"
    template = templateEnv.get_template(template_file)
    output = template.render(
                script=script, 
                div=div,
                protein_id=request_prot,
                protein_desc='pippo',
                start_hl=5,
                start_r=3
                )
    open('test3.html','w').write(output)       

    