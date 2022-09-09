function cmap = twocolor(startcolor,endcolor,nbins)

rvals = linspace(startcolor(1),endcolor(1),nbins).';
gvals = linspace(startcolor(2),endcolor(2),nbins).';
bvals = linspace(startcolor(3),endcolor(3),nbins).';
cmap = [rvals, gvals, bvals];

end