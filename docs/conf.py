project = 'UH FEL Simulation'
author = 'Eremey Valetov'
copyright = '2026, Eremey Valetov'

extensions = ['myst_parser']
myst_enable_extensions = ['dollarmath', 'amsmath']
source_suffix = {'.md': 'markdown'}
master_doc = 'index'

html_theme = 'alabaster'
html_theme_options = {
    'description': 'Documentation for the UH MkV FEL simulation codes.',
    'fixed_sidebar': True,
}

exclude_patterns = ['_build']

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '11pt',
}

latex_documents = [
    (master_doc, 'uh-fel-simulation.tex', 'UH FEL Simulation',
     'Eremey Valetov', 'manual'),
]
