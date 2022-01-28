# rlc color_b.py version 8.0
# Copyright (c) 2004 Robert L. Campbell
# added user defined colors for 3 color ramp- Mark A. Wall
# added selection string to color names to avoid overlap between different selections

import colorsys,sys,re
from pymol import cmd

# main function called from within PyMOL
def color_b(selection='all', item='b', mode='hist', gradient='bgr', nbins=11, sat=1.0, value=1.0, minimum='', maximum='', user_rgb='', debug=0):
  """

  AUTHOR

    Robert L. Campbell with enhancements from James Stroud

  USAGE

    color_b selection='sel',item='b', 'q', 'partial_charge' or 'formal_charge'
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' or
      'rw' or 'wr' or 'ry' or 'yr' or 'gw' or 'wg' or 'bw' or wb' or
      'gy' or 'yg' or 'gray' or 'reversegray' or 'user'

      mode='hist' or 'ramp' (default is 'hist')

      [minimum=''],[maximum=20.],

      nbins=11, sat=1.0, value=1.0,

      user_rgb = '(r1,g1,b1,r2,g2,b2,r3,g3,b3') [for use with gradient=user]

      The "item" argument allows specifying 'b', 'q', 'index', 'partial_charge'
      or 'formal_charge'as the item to color on.  The "color_q" function
      is really just the same as "color_b item=q".  Using item=index is
      similar to using the built-in "spectrum" command.

      This function allows coloring of a selection as a function of
      B-value or occupancy, following a gradient of colours.  The
      gradients can be:

      'bgr': blue -> green   -> red
      'rgb': red  -> green   -> blue
      'bwr': blue -> white   -> red
      'rwb': red  -> white   -> blue
      'bmr': blue -> magenta -> red
      'rmb': red  -> magenta -> blue
      'rw' : red -> white
      'wr' : white -> red
      'ry' : red -> yellow
      'yr' : yellow -> red
      'gw' : green -> white
      'wg' : white -> green
      'bw' : blue -> white
      'wb' : white -> blue
      'gy' : green -> yellow
      'yg' : yellow -> green
      'gray' : black -> white
      'reversegray' : white -> black
      'user' : user defined in this script

      ('rainbow' and 'reverserainbow' can be used as synonyms for
      'bgr' and 'rgb' respectively and 'grey' can be used as a synonym for 'gray').

      User-defined gradients are entered on the command line in
      parentheses as either integers between 0 and 255 or floats between
      0 and 1.  If any one value is larger than 1, then it is assumed
      that all are being entered as integers between 0 and 255.  Hence one can type:

      color_b selection, gradient=user, user_rgb=(0,0,1, 0,.5,1., 1.,.5,0.)

        or

      color_b selection, gradient=user, user_rgb=(0,0,255, 0,128,255, 255,128,0.)

      The division of B-value ranges can be in either of two modes: 'hist' or
      'ramp'. 'hist' is like a histogram (equal-sized B-value increments
      leading to unequal numbers of atoms in each bin). 'ramp' as a ramp
      of B-value ranges with the ranges chosen to provide an equal number
      of atoms in each group.

      You can also specify the lower or upper limits of the data used to determine
      the color bins (minimum,maximum). e.g. color_b my_molecule, minimum=15., maximum=25.

      You can also specify the saturation and value (i.e. the "s" and "v"
      in the "HSV" color scheme) to be used for the gradient. The defaults
      are 1.0 for both "sat" and "value".

      In the case of the gray scale gradients, "sat" sets the minimum intensity
      (normally black) and "value" sets the maximum (normally white)

    usage:
      from within PyMOL do "run color_b.py" to load the function definition.
      Then you can use for example:

          color_b (c. a | c. b),mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1.

      to color chains A and B with the Blue-White-Red gradient in 30 colors of equal
      numbers of atoms in each color.
  """

  nbins=int(nbins)
  sat=float(sat)
  value=float(value)
# make sure sat and value are in the range 0-1.0
  sat = min(sat, 1.0)
  sat = max(sat, 0.0)
  value = min(value, 1.0)
  value = max(value, 0.0)
  debug = int(debug)
  if gradient == 'user' and user_rgb == '':
    user_rgb = '50,50,195, 245,245,20, 255,20,20'

# make sure lowercase
  gradient.lower()
  mode.lower()

# Sanity checking
  if nbins == 1:
    print("\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n")
    nbins=11

  if mode not in ('hist','ramp'):
    print("\n     WARNING: Unknown mode ",mode, "    ----->   Nothing done.\n")
    return
  elif gradient not in ('bgr','rgb','rainbow','reverserainbow','bwr','rwb','user',
                        'bmr','rmb','rw','wr','ry','yr','gw','wg','bw','wb','gy',
                        'yg','gray','grey','reversegray','reversegrey','user'):
    print("\n     WARNING: Unknown gradient: ",gradient, "    ----->   Nothing done.\n")
    return

  #print("MODE, GRADIENT, NBINS:", mode,gradient, nbins)

# get list of B-factors from selection

# contents of "m.atom[i]": 'b', 'chain', 'coord', 'defaults',
#'elec_radius', 'flags', 'formal_charge', 'get_implicit_valence',
#'get_mass', 'get_number', 'get_signature', 'has', 'hetatm', 'id',
#'in_same_residue', 'index', 'name', 'new_in_residue', 'numeric_type',
#'partial_charge', 'q', 'resi', 'resi_number', 'resn', 'segi', 'ss',
#'stereo', 'symbol', 'u_aniso', 'vdw']

  m = cmd.get_model(selection)
  sel = []
  b_list = []

  if len(m.atom) == 0:
    print("Sorry, no atoms selected")

  else:
    if item == 'b':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].b)
    elif item == 'q':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].q)

    elif item == 'partial_charge':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].partial_charge)

    elif item == 'formal_charge':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].formal_charge)

    elif item == 'index':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].index)

    elif item == 'resi':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].resi_number)

    else:
      print("Not configured to work on item %s" % item)
      return

    max_b = float(max(b_list))
    min_b = float(min(b_list))
    #print("Minimum and Maximum B-values: ", min_b, max_b)

    if mode == 'hist':

# check if minimum or maximum was specified and use the entered values
      if minimum != '':
        min_b = float(minimum)
      if maximum != '':
        max_b = float(maximum)
      # histogram:
      # color in bins of equal B-value ranges
      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
      bin_width = (max_b - min_b)/nbins*1.0
      if item == 'b':
        b_count = get_signle_atom_props(selection)
        b_factors = set(list(b_count['bfactors']))
        for b_value in b_factors:
            temp = selection + " and (%s  = %4.4g" % (item,b_value) + ")"
            sel.append(temp)
        
# call the function to create the gradient which returns a list of colours
    colours = make_gradient(sel,gradient,len(sel),sat,value, user_rgb,debug)

# do the colouring now
    for j in range(len(sel)):
      #print("Color select: ",sel[j])
      cmd.color(colours[j],sel[j])
  #    print(j,colours[j],sel[j])

def get_signle_atom_props(selection):
    stored_b = []
    myspace = {'bfactors': []}
    cmd.iterate('(all)', 'bfactors.append(b)', space=myspace)
    return myspace

# function for creating the gradient
def make_gradient(sel,gradient,nbins,sat,value,user_rgb,debug=0):
  if gradient == 'bgr' or gradient == 'rainbow':
    gradient = 'bgr'
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique
      # for the selection
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # coldesc.append('col' + str(sel[j]) + str(j))

      # create colors using hsv scale (fractional) starting at blue(.6666667)
      # through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1"
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      #cmd.set_color("col" + gradient + str(j),col[j])
      #if debug:
        #print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])
      #print(col[j],"defined as ", "col"+str(j))

#  if debug:
#    for j in range(nbins):
#      print("colour #:",j,"colour RGB triplet: ",col[j])

  #print(coldesc)
# return the gradient as a list of colors named by their gradient & index (i.e. colbw0,colbw1,colbw2,...)
  return coldesc

# allow calling without parentheses: color_hist_b [selection=], [mode= ],[gradient= ],[nbins= ]
cmd.extend("color_b",color_b)
