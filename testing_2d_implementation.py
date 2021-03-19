import time
import region_growing_python as rgp

seed = (11,45,35)

start = time.time()
seg = rgp.grow(vol, seed, 5)
stop = time.time()

print("Elapsed time: %.3f seconds." % (stop - start))
print("Errors: %d" % np.logical_xor(w <= -11, seg).sum())

src = mlab.pipeline.scalar_field(seg.astype(np.float))
mlab.pipeline.iso_surface(src, contours=[0.5], opacity=0.5)
mlab.show()
