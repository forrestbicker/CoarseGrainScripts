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

    measurement_blueprint = list(measurement_dict.keys())[0]
    measurement_type = MesType(len(measurement_blueprint))
    measurement_name = '-'.join(measurement_blueprint)

    data = measurement_dict[measurement_blueprint]
    aggregate_values = []
    for entry in data.values():
        aggregate_values += entry['values']
    bin_size = 0.01
    hist = Histogram(measurement_type, aggregate_values, bin_size, name=measurement_name)
    hist.clear_empty_bins()

    x_data = hist.get_floors()
    y_data = hist.get_boltzes()

    original_fig, k, x0, c = generate_figure(x_data, y_data, hist.name)

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
                    figure=original_fig,
                ),
            ])
        ]
    )


    @app.callback(
        [
            Output('fitting-graph', 'figure'),
            Output('results', 'children'),
            Output('page-number', 'children'),
        ],
        [
            Input('bin-slider', 'value'),
            Input('next-button', 'n_clicks'),
            Input('back-button', 'n_clicks'),
            Input('fitting-graph', 'selectedData'),
        ],
        prevent_initial_call=True,
    )
    def display_selected_data(bin_size, _, _2, selectedData):

        global measurement_number

        ctx = dash.callback_context
        triggered_input = ctx.triggered[0]['prop_id'].split('.')[0]

        if triggered_input == 'next-button':
            if measurement_number < len(measurement_dict) - 1:
                measurement_number += 1

        if triggered_input == 'back-button':
            if measurement_number > 0:
                measurement_number -= 1

        measurement_blueprint = list(measurement_dict.keys())[measurement_number]
        measurement_type = MesType(len(measurement_blueprint))
        measurement_name = '-'.join(measurement_blueprint)

        #  when need to redraw data from the dictionary
        if triggered_input in ['next-button', 'back-button', 'bin-slider'] or selectedData is None:
            data = measurement_dict[measurement_blueprint]
            aggregate_values = []
            for entry in data.values():
                aggregate_values += entry['values']

            hist = Histogram(measurement_type, aggregate_values, bin_size, name=measurement_name)
            hist.clear_empty_bins()
            x_data = hist.get_floors()
            y_data = hist.get_boltzes()
            biggest_bin = hist.get_biggest(1)[0]
            vertex = (biggest_bin.floor, biggest_bin.boltz())
        # otherwise, draw data from the currently selected data
        # TODO: nail down exactly how the selected_data variable works,
        # will it default to all data or null?
        # what format will it be in?
        else:
            x_data = [point['x'] for point in selectedData['points']]
            y_data = [point['y'] for point in selectedData['points']]
            biggest_bin_ix = y_data.index(min(y_data))
            vertex = (x_data[biggest_bin_ix], y_data[biggest_bin_ix])

        fig, k, x0, c = generate_figure(x_data, y_data, measurement_name, vertex=vertex)
        equation = f'Current line of fit: y = {k: .3f}(x - {x0: .3f}) + {c: .3f}'
        progress_str = f'Currently viewing {measurement_type.name} {measurement_name}. Measurement {measurement_number + 1} / {len(measurement_dict)}'

        potentials[measurement_blueprint] = (k, x0, c)

        return fig, equation, progress_str


