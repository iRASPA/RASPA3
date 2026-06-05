import ipywidgets as widgets
from IPython.display import display
import py3Dmol
import io
from ase.io import read
from ase.spacegroup import crystal
from ase.io import write


def MFI(x,y,z):
  number_of_unit_cells = (x, y, z)
  a, b, c = 20.022, 19.899, 13.383
  alpha, beta, gamma = 90, 90, 90
  
  species = ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si", \
             "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O"]
  coords = [
      [ 0.42238,  0.0565,  -0.33598], [ 0.30716,  0.02772, -0.1893 ], [ 0.27911,  0.06127,  0.0312 ], [ 0.12215,  0.06298,  0.0267 ], [ 0.07128,  0.02722, -0.18551], [ 0.18641,  0.05896, -0.32818],
      [ 0.42265, -0.1725,  -0.32718], [ 0.30778, -0.13016, -0.18548], [ 0.27554, -0.17279,  0.03109], [ 0.12058, -0.1731,   0.02979], [ 0.07044, -0.13037, -0.182  ], [ 0.18706, -0.17327, -0.31933],
      [ 0.3726,   0.0534,  -0.2442 ], [ 0.3084,   0.0587,  -0.0789 ], [ 0.2007,   0.0592,   0.0289 ], [ 0.0969,   0.0611,  -0.0856 ], [ 0.1149,   0.0541,  -0.2763 ], [ 0.2435,   0.0553,  -0.246  ],
      [ 0.3742,  -0.1561,  -0.2372 ], [ 0.3085,  -0.1552,  -0.0728 ], [ 0.198,   -0.1554,   0.0288 ], [ 0.091,   -0.1614,  -0.0777 ], [ 0.1169,  -0.1578,  -0.2694 ], [ 0.2448,  -0.1594,  -0.2422 ],
      [ 0.3047,  -0.051,   -0.1866 ], [ 0.0768,  -0.0519,  -0.1769 ], [ 0.4161,   0.1276,  -0.3896 ], [ 0.4086,  -0.0017,  -0.4136 ], [ 0.402,   -0.1314,  -0.4239 ], [ 0.1886,   0.1298,  -0.3836 ],
      [ 0.194,    0.0007,  -0.4082 ], [ 0.1951,  -0.1291,  -0.419  ], [-0.0037,   0.0502,  -0.208  ], [-0.004,   -0.1528,  -0.2078 ], [ 0.4192,  -0.25,    -0.354  ], [ 0.1884,  -0.25,    -0.3538 ],
      [ 0.2883,  -0.25,     0.0579 ], [ 0.1085,  -0.25,     0.0611 ]
  ]
  
  structure = crystal(symbols=species, basis=coords, spacegroup=62, cellpar=[a, b, c, alpha, beta, gamma]) * number_of_unit_cells
  
  cif_io = io.BytesIO()
  write(cif_io, structure, format='cif')
  cif_str = cif_io.getvalue().decode('utf-8')

  return cif_str


def ITQ_29(x,y,z):
  number_of_unit_cells = (x, y, z)
  a, b, c = 11.8671, 11.8671, 11.8671
  alpha, beta, gamma = 90, 90, 90
  
  species = ["Si", "O", "O", "O"]
  coords = [
      [0.3683, 0.1847, 0],
      [0.5,    0.2179, 0],
      [0.2939, 0.2939, 0],
      [0.3429, 0.1098, 0.1098]
  ]
  
  structure = crystal(species, basis=coords, spacegroup=221, cellpar=[a, b, c, alpha, beta, gamma]) * number_of_unit_cells

  cif_io = io.BytesIO()
  write(cif_io, structure, format='cif')
  cif_str = cif_io.getvalue().decode('utf-8')

  return cif_str

def create_molecular_movie(view, pdb_path, framework=None, interval=200):
  def pos(arr):
      return {'x': float(arr[0]), 'y': float(arr[1]), 'z': float(arr[2])}

  box = read(pdb_path).get_cell()
  a, b, c = box[0], box[1], box[2]
  v0 = [0, 0, 0]; v1, v2, v3 = a, b, c; v12, v13, v23 = a + b, a + c, b + c; v123 = a + b + c
  edges = [(v0, v1), (v0, v2), (v0, v3), (v1, v12), (v1, v13), (v2, v12), (v2, v23), (v3, v13), (v3, v23), (v12, v123), (v13, v123), (v23, v123) ]
  for start_node, end_node in edges:
      view.addCylinder({'start': pos(start_node), 'end': pos(end_node), 'radius': 0.2, 'color': 'white', 'fromCap': 1, 'toCap': 1})
  
  pdb_string = open(pdb_path, 'r').read()
  view.addModelsAsFrames(pdb_string,'multimodelpdb')

  if framework:
    view.addModel(framework, 'cif')

  view.setStyle({}, {})
  view.addStyle({'model': 0}, {'sphere': {'scale': 1.0,}})
  if framework:
    view.setStyle({'model': 1}, {'stick': {'radius': 0.15}, 'sphere': {'scale': 0.15}})

  view.setCameraParameters({'orthographic': True})

  initial_orientation = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]

  play_btn = widgets.Button(description="Play", icon="play")
  stop_btn = widgets.Button(description="Stop", icon="stop")
  reset_btn = widgets.Button(description="Reset", icon="history")
  
  def play_movie(b):
      view.animate({'loop': 'once', 'interval': interval, 'reps': 1})
      view.update()
  
  def stop_movie(b):
      view.stopAnimate()
      view.update()

  def reset_view(b):
    view.stopAnimate()
    view.setView(initial_orientation)
    view.setFrame(0)
    view.zoomTo()
    view.update()
  
  play_btn.on_click(play_movie)
  stop_btn.on_click(stop_movie)
  reset_btn.on_click(reset_view)

  display(widgets.HBox([play_btn, stop_btn, reset_btn]))
  view.setView(initial_orientation)
  view.zoomTo()
