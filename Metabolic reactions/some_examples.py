# Importing the cobra library
import cobra
from cobra import Model, Reaction, Metabolite

# step 1 : Constructing a metabolic model

model = Model('example_model')

# set of reactions
reactions = dict()
for i in range(1,8):
    reactions["R"+str(i)] = Reaction("R"+str(i))
    reactions["R"+str(i)].name = "R"+str(i)
    reactions["R"+str(i)].subsystem = 'Body'
    reactions["R"+str(i)].lower_bound = 0.  # This is the default
    reactions["R"+str(i)].upper_bound = 1000.  # This is the default
    #reactions["R"+str(i)].objective_coefficient = 0. # this is the default
    print ("Creating R"+str(i))
# set of metabolites
A = Metabolite('A', name='A', compartment='c')
B = Metabolite('B', name='B', compartment='c')
C = Metabolite('C', name='C', compartment='c')
D = Metabolite('D', name='D', compartment='c')
E = Metabolite('E', name='E', compartment='c')
F = Metabolite('F', name='F', compartment='c')

reactions["R1"].add_metabolites({A:-2, B:-1, C: 1})
reactions["R2"].add_metabolites({B:-1, C:-1, D: 1, E: 2})
reactions["R3"].add_metabolites({A:-4, D: 1})
reactions["R4"].add_metabolites({A:1})
reactions["R5"].add_metabolites({B:1})
reactions["R6"].add_metabolites({D: -1})
reactions["R7"].add_metabolites({E: -1})


reactions["R4"].lower_bound = 1
reactions["R4"].upper_bound = 6
reactions["R7"].lower_bound = 1
reactions["R7"].upper_bound = 4

for i in range(1,8):
    model.add_reaction(reactions["R"+str(i)])


# setting the objective

model.reactions.R4.objective_coefficient=1
model.reactions.R7.objective_coefficient=1

# or alternative solution
model.objective= model.reactions.R4.flux_expression + model.reactions.R7.flux_expression

# Display some informations on the model (part 1)
print ("Number of metabolites = ",len(model.metabolites))
print ("Number of reactions = ",len(model.reactions))

# Display some informations on the model (part 2)

print ("Reaction 2 is ", model.reactions[2].build_reaction_string(), "lb = ", model.reactions[2].lower_bound, "ub = ", model.reactions[2].upper_bound)

# Question : Write a loop in order to get all the reactions

# Notice than we can iterate on reactions (or metabolites or genes,...)

for r in model.reactions:
    print(r.name, r.build_reaction_string())



# Perform a FBA optimization

# performing a FBA is done by using the optimize() function that
# construct the appropriate LP problem and call the solver.

f=model.optimize()
print ("optimal value=",f.objective_value)
print ("corresponding fluxes = ",f.fluxes)

# a nicer view with pandas library
import pandas
pandas.DataFrame.from_dict({"fluxes": f.fluxes})

# Question: Get the minimal and maximal values of all the fluxes in the solution space

# for that, we need to change the objective and perform FBAs

for r in model.reactions:
    # for the maximal value of flux v_r (r is the reaction), we set obj(v)=v_r like this
    model.objective = r.flux_expression
    # the we optimize
    vr_max = model.optimize()
    # for the minimal value of flux v_r (r is the reaction), we set obj(v)=-v_r
    # and use the fact that min(v_r) = -max(-v_r)
    model.objective = -r.flux_expression
    # the we optimize
    vr_min = model.optimize()
    print(r.id, "min : ", -vr_min.objective_value," , max : ",vr_max.objective_value)


# Extra ! A plot of the solution space by picking randomly some fluxes
# belonging to the solution space


nbpoints=1000
PointsR4 = list()
PointsR7 = list()
PointsR1 = list()
PointsR2 = list()
fluxes_sample = cobra.sampling.sample(model,nbpoints)

# We are only interested by V_4 and V_7
for i in range(nbpoints):
    PointsR4.append(fluxes_sample.get("R4")[i])
    PointsR7.append(fluxes_sample.get("R7")[i])
    PointsR1.append(fluxes_sample.get("R1")[i])
    PointsR2.append(fluxes_sample.get("R2")[i])


# the following plot will be an approximation of the solution space

import matplotlib.pyplot as plt
%matplotlib inline

plt.plot(PointsR4,PointsR7,"ro")
plt.show()

plt.plot(PointsR1,PointsR4,"bo")
plt.show()

plt.plot(PointsR1,PointsR2,"go")
plt.show()


import pandas
import cobra.test

# You can use this model from the tutorial
#model = cobra.test.create_test_model("ecoli")

# or this one that must be in the home folder
model = cobra.io.read_sbml_model("e_coli_core.xml")


# Question: Execute a for loop on all the reactions and print their id and objective coefficient.

for r in model.reactions:
    if not r.objective_coefficient == 0:
        print(r.id, r.objective_coefficient)

# notice that this verification is important because several models may not contains an objective function
# this test must be done before trying to compute a FBA of course

# Question: how many reactions have a zero 0 for an optimal solution

# 1. perform a FBA

fba = model.optimize()

# 2. then count the number of zeros in fba.fluxes
zerofluxes = 0
for r in model.reactions:
    vr = fba.fluxes.get(r.id)
    if vr == 0:
        zerofluxes = zerofluxes + 1
print("Reaction with a null flux : ",zerofluxes)

# Perform a FVA study. Classify the reactions into 3 classes.
# Blocked reactions : ----------[0]----------
# Essential reactions : ----------0--[------]-- or --[------]--0----------
# Alternative reactions : --------[--0---]-------

# 1. Computes a FVA (depending on alpha)

alpha = float(input("Enter the fraction of optimum value : "))
fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=alpha)

# 1. Computes the blocked and essential reactions
blocked = 0
essential = 0
for r in model.reactions:
    vr_min = fva.get("minimum").get(r.id)
    vr_max = fva.get("maximum").get(r.id)
    if (abs(vr_min)<0.00000000001 and abs(vr_max)<0.00000000001):
        blocked = blocked + 1
    if vr_min>0.00000000001 or vr_max<-0.00000000001:
        essential = essential + 1

alternative = len(model.reactions)-blocked-essential

print(blocked, "blocked reactions")
print(essential, "essential reactions")
print(alternative, "alternative reactions")

data = list()

for i in range(1000):
    model.reactions.EX_glc__D_e.lower_bound = -i
    model.reactions.EX_glc__D_e.upper_bound = 0
    fba = model.optimize()
    data.append(fba.objective_value)

import matplotlib.pyplot as plt
plt.plot(data)
plt.show()


model = cobra.io.read_sbml_model("e_coli_core.xml")
fba_ref = model.optimize()

for i in model.reactions:
    tmplb = i.lower_bound
    tmpub = i.upper_bound
    i.lower_bound = 0
    i.upper_bound = 0
    fba = model.optimize()
    if fba.objective_value < 0.9*fba_ref.objective_value:
        print(i.id, fba.objective_value)
    i.lower_bound = tmplb
    i.upper_bound = tmpub

for r in model.reactions:
    print(r.build_reaction_string(), r.lower_bound, r.upper_bound)
