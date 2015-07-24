#ifndef __markers_h
#define __markers_h

PetscErrorCode allocateMarkers(Problem *);
PetscErrorCode destroyMarkers(MarkerSet *,Options *);
void resetMarkers( MarkerSet *,Options *);
//PetscErrorCode distributeMarkers( Markers *, PetscInt, PetscInt, PetscInt , PetscInt, PetscScalar, PetscScalar);
PetscErrorCode distributeMarkersUniformInDomain( MarkerSet *, Options *, GridData *);
PetscErrorCode distributeMarkersUniformInCells( MarkerSet *, Options *, GridData *);
PetscErrorCode projectMarkersNodesAll(MarkerSet *, GridData *, NodalFields *, Materials * ,Options *);
PetscErrorCode projectMarkersNodesFromScalar(MarkerSet *, GridData *, PetscScalar *, Vec);
PetscErrorCode projectMarkersCellCentersFromScalar(MarkerSet *, GridData *, PetscScalar *, Vec);
PetscErrorCode projectVelocityToMarkers(MarkerSet *, GridData *, NodalFields * );
PetscErrorCode checkMarkersOutsideDomain( MarkerSet *, PetscScalar , PetscScalar );
PetscErrorCode checkPlasticYielding(GridData*, MarkerSet *, Materials *, PetscScalar,PetscInt *,Options *, PetscScalar);
void findMarkerCells(MarkerSet *, GridData *);
void printMarkersAllASCIIMatlab( MarkerSet*, PetscInt);
PetscErrorCode exchangeMarkers( MarkerSet* , GridData *,Options *);
PetscErrorCode checkMarkerDensity( Problem *, MarkerSet *, GridData *, Options *,PetscRandom );
PetscErrorCode advectMarkers(MarkerSet *, GridData *, PetscScalar );
PetscErrorCode advectMarkersRK(MarkerSet *, NodalFields *, GridData *, Options *, PetscScalar );
PetscErrorCode projectNodalFieldToMarkersS(NodalFields *, Vec, MarkerSet *, PetscScalar *, GridData *);
void resetMarker( Marker * );

#endif
