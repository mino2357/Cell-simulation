border a(t=0, 3.0){x=-1.5+t; y=-1.5;};
border b(t=0, 3.0){x=1.5; y=-1.5+t;};
border c(t=0, 3.0){x=1.5-t; y=1.5;};
border d(t=0, 3.0){x=-1.5; y=1.5-t;};
int cut = 8;
mesh Th = buildmesh(a(cut)+b(cut)+c(cut)+d(cut));
savemesh(Th, "TT.msh");
