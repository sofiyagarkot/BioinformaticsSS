from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='histamineR2.pir',
              knowns='3rze', sequence='HHR2',
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
