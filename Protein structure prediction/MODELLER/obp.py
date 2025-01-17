from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='obp.pir',
              knowns='4run', sequence='Q9NPH6',
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
