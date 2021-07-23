import plotly.graph_objects as go
from ipywidgets import widgets
from IPython.display import display
import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import pandas as pd
from dash.dependencies import Input, Output
from dash_extensions import Download


import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd

from BINAnalysis import Histogram
from util import MesType, generate_figure

app = dash.Dash(__name__)

def wrap(component):
    return html.Div(
        children=[
            html.Div(
                children=component,
                className='card'
            )
        ],
        className='wrapper'
    )

measurement_number = 0
potentials = {}

def startup_manual_refining(measurement_dict, u):


    app.layout = html.Div(
        children=[
            wrap([
                html.H1('Curve Fitting Tool'),
                html.P('Slide the slider to adjust the size of the bins that the raw measurements are clustered into. Select a region of interest on the graph to adjust the subset of data to which the curve is fitted. When you are happy with the results, click NEXT to view the next graph or EXPORT to write your parametization data to a file.'),
                html.H2('Adjust Bin Size'),
                dcc.Slider(
                    id='bin-slider',
                    min=0.001,
                    max=1,
                    step=0.001,
                    value=0.01,
                ),
                html.H2(id='page-number'),
                html.Button(id='back-button', n_clicks=0, children='Back'),
                html.Button(id='next-button', n_clicks=0, children='Next'),
                html.Button(id='export-button-1', n_clicks=0, children='Export Params'),
                Download(id='export-download-1'),
                html.Button(id='export-button-2', n_clicks=0, children='Export .top'),
                Download(id='export-download-2'),
                html.H2(id='results'),
                dcc.Graph(
                    id='fitting-graph',
                ),
            ])
        ]
    )




