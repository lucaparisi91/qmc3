# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from dash.dependencies import Input, Output
from jupyter_dash import JupyterDash

def plotConfiguration(data):
        data["particle"]=data["particle"].astype(str)
        graphs=[]
        for particle,particleData in data.groupby("particle"):
            graphs.append(go.Scatter3d( x=particleData['x'], y=particleData['y'], z=particleData['z'] ,name=str(particle)))
        fig=go.Figure(data=graphs)
        return fig
        

class configurationsApp:
    
    def __init__(self,filename):
        external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

        self.app = JupyterDash(__name__, external_stylesheets=external_stylesheets)
        fig=go.Figure(data=[])

        try:
            conf=pd.read_csv(filename,sep=" ",index_col=False)
            fig=plotConfiguration(conf)

        except FileNotFoundError:
            print ("File not found")

        self.app.layout = html.Div(children=[
            html.H1(children='Visualization of PIMC simulation'),

    
            html.H2(children='''
                Configurations
            '''),

            dcc.Graph(
                id='last-configuration',
                figure=fig
            )
            ])