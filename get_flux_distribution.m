%this code requires the Cobra Toolbox to function properly
%First, a branched pathway is constructed.
%Second, a possible flux distribution is calculated through this pathway
initCobraToolbox

%pre-allocate
model.rxns      =       {};
model.mets      =       {};
model.S         =       [];
model.rev       =       [];
model.lb        =       [];
model.ub        =       [];
model.c         =       [];
model.rules     =       {};
model.genes     =       {};

%setup model reactions
model           =       addReaction(model,'r1',' -> orn');
model           =       addReaction(model,'r2',' -> cbp');
model           =       addReaction(model,'r3',' orn + cbp -> arg');
model           =       addReaction(model,'r4',' cbp -> utp');
model           =       addReaction(model,'r5',' arg -> ');
model           =       addReaction(model,'r6',' utp -> ');

%define rxn bounds and solve system
model           =       changeRxnBounds(model,'r5',0.958,'b');
model           =       changeRxnBounds(model,'r6',0.467,'b');
model           =       changeObjective(model,'r3');
solution        =       optimizeCbModel(model,'max');

%get steady state flux
flux_ss         =       solution.x;

%check flux for steady state
SS              =       model.S*flux_ss;