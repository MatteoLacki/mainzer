%load_ext autoreload
%autoreload 2
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input, State
import dash_table

import pandas as pd
import io
import base64
import datetime

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
# app = dash.Dash(__name__)

styles = {
    "upload": {
            'width': '30%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        }
}

def ask_for(var_name, message, type="text", default=None):
    return html.Div([html.Label(message),
                     dcc.Input(id=var_name, type="number", value=default)])


app.layout = html.Div([
    html.H1("Welcome to Lipido!"),
    html.H2("we like it top down"),
    dcc.Upload(
        id='molecules_upload',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select a CSV/Excel file describing Molecules')
        ]),
        style=styles["upload"],
        multiple=True
    ),
# TODO NEED TO CONSIDER HOW TO UPLOAD A BINARY FILE LIKE A ZIP

    html.Div(id='output-data-upload',
             style={'width':"50%"}),
    dcc.Upload(
        id='spectrum-upload',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select a Spectrum File [csv|mzml|mzxml]')
        ]),
        style=styles["upload"],
        multiple=True
    ),
    html.Div(id='output-spectrum-upload',
             style={'width':"50%"}),

    html.Div(id='spectrum_upload',
             style={'width':"50%"}),
    ask_for("min_intensity", "Select Minimal Intensity: ", "number", 0),
    ask_for("max_lipid_mers", "Maximal Number of Lipid Mers: ", "number", 4),
    ask_for("min_protein_charge", "Minimal Protein Chargestate: ", "number", 1),
    ask_for("max_protein_charge", "Maximal Protein Chargestate: ", "number", 10),
    ask_for("min_lipid_charge", "Minimal Charge on Free Lipid Cluster: ", "number", 1),
    ask_for("max_lipid_charge", "Maximal Charge on Free Lipid Cluster: ", "number", 1),
    ask_for("isotopic_coverage", "IsoSpec probability coverage [0<x<1]: ", "number", .95),
    ask_for("isotopic_bin_size", "IsoSpec bin size in Thomsons [0<x<1]: ", "number", .1),
    ask_for("neighbourhood_thr", "Neighbourhood buffer size in Thomsons [0<x]: ", "number", 1.1),
    ask_for("underfitting_quantile", "Single molecule underfit quantile [0<x]: ", "number", 0.5),
    dcc.Checklist(
        options=[
            {'label': 'Deconvolve?', 'value': 'deconvolve'},
        ],
        value=['NYC', 'MTL']
    ),
    dcc.Textarea(
        id='comment',
        value='Maybe a sample description?',
        style={'width': '30%', 'height': 100},
    ),
    html.Div([
        html.Button('Run Lipido!', id='submit', n_clicks=0),
    ]),
    html.Div(id='container-button-basic',
             children='Enter a value and press submit')
])


def parse_upload(contents, filename, date):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if '.csv' in filename:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            df = pd.read_excel(io.BytesIO(decoded))
        elif 'mzML' in filename:
            import pyteomics.mzml
            import pyteomics.mzxml
            test = pyteomics.mzml.read(decoded.decode('utf-8'))
            data = next(test)
            df = pd.DataFrame({"mz":data['m/z array'],
                               "intensity":data['intensity array']})
            print(df)
            # data['m/z array'], data['intensity array']
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return dash_table.DataTable(data=df.to_dict('records'),
                                columns=[{'name': i, 'id': i} for i in df.columns],
                                style_cell={
                                    'overflow': 'hidden',
                                    'textOverflow': 'ellipsis',
                                    'maxWidth': 0,
                                },
                                tooltip_data=[
                                    {
                                        column: {'value': str(value), 'type': 'markdown'}
                                        for column, value in row.items()
                                    } for row in df.to_dict('records')
                                ],
                                tooltip_duration=None)



@app.callback(Output('output-data-upload', 'children'),
              Input('molecules_upload', 'contents'),
              State('molecules_upload', 'filename'),
              State('molecules_upload', 'last_modified'))
def update_output(list_of_contents, list_of_names, list_of_dates):
    print(list_of_contents, list_of_names, list_of_dates)
    if list_of_contents is not None:
        children = [
            parse_upload(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children

@app.callback(Output('output-spectrum-upload', 'children'),
              Input('spectrum-upload', 'contents'),
              State('spectrum-upload', 'filename'),
              State('spectrum-upload', 'last_modified'))
def update_output2(list_of_contents, list_of_names, list_of_dates):
    print(list_of_contents, list_of_names, list_of_dates)
    if list_of_contents is not None:
        children = [
            parse_upload(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children


@app.callback(
    dash.dependencies.Output('container-button-basic', 'children'),
    [
        Input("submit", "n_clicks")
    ],
    [
        # State("molecules_csv", "value"),
        # State("spectrum", "value"),
        State("min_intensity", "value"),
        State("max_lipid_mers", "value"),
        State("min_protein_charge", "value"),
        State("max_protein_charge", "value"),
        State("min_lipid_charge", "value"),
        State("max_lipid_charge", "value"),
        State("isotopic_coverage", "value"),
        State("isotopic_bin_size", "value"),
        State("neighbourhood_thr", "value"),
        State("underfitting_quantile", "value")
    ]
)
def number_render(n_clicks,
                  # molecules_csv,
                  # spectrum,
                  min_intensity, 
                  max_lipid_mers,
                  min_protein_charge,
                  max_protein_charge,
                  min_lipid_charge,
                  max_lipid_charge,
                  isotopic_coverage,
                  isotopic_bin_size,
                  neighbourhood_thr,
                  underfitting_quantile):
    print(locals())
    return locals()



if __name__ == '__main__':
    app.run_server()