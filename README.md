# UNITN-CS-2020-SignalVideoImaging-project

http://notmatthancock.github.io/2017/10/09/region-growing-wrapping-c.html

## Pseudocode
```
method grow(image, seed, neighborhood_size)

Initialize 'segmentation' to boolean volume, same shape as 'image'.
Initialize 'checked' to a boolean volume, same shape as 'segmentation'.
Initialize empty stack, 'needs_check'.

Set 'segmentation' and 'checked' at 'seed' to true.
Add neighbor coordinates of 'seed' to 'needs_check'.

while 'needs_check' is not empty:
    Pop 'point' from 'needs_check'.
    Set 'checked' at 'point' to true.

    if the average of 'image' over 'neighborhood_size' distance from 'point'
       is greater than or equal to the value of 'image' at 'point'
    then
        Set 'segmentation' to true at 'point'.
        Add neighbor coordinates of 'point' to 'needs_check'.
    end if
end while

return 'segmentation'
end method grow
```

