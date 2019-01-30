m = document.getModel();
printAnnotation(m);

for i in range(0, m.getNumReactions()):
  re = m.getReaction(i);
  printAnnotation(re);

    # SpeciesReference (Reacatant)
  for j in range(0, re.getNumReactants()):
    rt = re.getReactant(j);
    if (rt.isSetAnnotation()):
      print("     ");
    printAnnotation(rt, rt.getSpecies());

  # SpeciesReference (Product) 
  for j in range(0, re.getNumProducts()):
    rt = re.getProduct(j);
    if (rt.isSetAnnotation()):
      print("     ");
    printAnnotation(rt, rt.getSpecies());

  # ModifierSpeciesReference (Modifiers)
  for j in range(0, re.getNumModifiers()):
      md = re.getModifier(j);
      if (md.isSetAnnotation()):
        print("     ");
      printAnnotation(md, md.getSpecies());

  # KineticLaw 
  if (re.isSetKineticLaw()):
    kl = re.getKineticLaw();
    if (kl.isSetAnnotation()):
      print("   ");
    printAnnotation(kl);

    # Parameter   
    for j in range(0, kl.getNumParameters()):
      pa = kl.getParameter(j);
      if (pa.isSetAnnotation()):
        print("      ");
      printAnnotation(pa);
  
# Species 
for i in range(0, m.getNumSpecies()):
  sp = m.getSpecies(i);
  printAnnotation(sp);
  
# Compartments 
for i in range(0, m.getNumCompartments()):
  sp = m.getCompartment(i);
  printAnnotation(sp);
