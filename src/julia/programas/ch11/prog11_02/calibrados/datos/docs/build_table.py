import yaml
import pandas as pd 


## Markdown table
with open("parameters_metadata.yml","r") as file:
    parametros = yaml.safe_load(file)

df = pd.DataFrame(
    [[v for k,v in param.items()] for param in parametros["parametros"]],
    columns = ["Parámetro", "Descripción", "Parámetro exógeno", "Parámetro Calibrado", "Target Calibración", "Fuente"]
)
print(df.to_markdown())

## Latex table

with open("parameters_metadata_latex.yml","r") as file:
    parametros = yaml.safe_load(file)

from jinja2 import Environment, FileSystemLoader

## Configuración de jinja 
file_loader = FileSystemLoader('.')
env = Environment(loader=file_loader)

template = env.get_template('template')


output = template.render(parametros = parametros["parametros"])

text_file = open("latex_parametros/dsolg_parametros.tex", "w")
text_file.write(output)
text_file.close()

