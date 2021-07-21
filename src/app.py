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
