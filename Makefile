default:
	(cd ann; make linux-g++)
	(cd OpenMesh; make)
	(cd pbrt/src; make)
	(cd registration/trimesh2; make)
	(cd registration/tps_alignment/src; make)
	(cd modeling; make)
	(cd sampler; make)
	(cd evaluator; make)

clean:
	(cd ann; make realclean)
	(cd OpenMesh; make clean)
	(cd pbrt/src; make clean)
	(cd registration/trimesh2; make spotless)
	(cd registration/tps_alignment/src; make clean)
	(cd modeling; make clean)
	(cd sampler; make clean)
	(cd evaluator; make clean)
