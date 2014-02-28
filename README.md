LIFbump
=======

leaky integrate and fire conductance-based (spiking) bump model

Model based on:

S. Vattikuti and C.C. Chow, 'A computational model for cerebral cortical dysfunction in Autism Spectrum Disorders', 
Biol Psychiatry 67:672-6798 (2010).

Please cite the above if you use this code for any published work. 

About:

This is a leaky integrate and fire bump model with simulated NMDA and GABA channels.
It is placed in the context of a cognitive switching activity in which some information is memorized (in working memory) 
at the beginning of the task and sometime later some new information (with some "distance" away from the inital datum) 
needs to be memorized. A bump is a confined region of elevated activity that persists after an external input has been extinguished. 
Such models have been used to emulate the dynamics of working memory based on single neuron recordings during delay, oculomotor tasks. 
The simulation shows how stability is affected by the balance of excitation (NMDA) and inhibition (GABA),
and by density. Density is altered by changing the neuron number but fixing the total field distance, such that an increase in 
neuron number leads to more densely packed neurons. Synaptic footprints and inputs are normalized. 
The network has periodic boundary conditions to avoid edge effects. 
