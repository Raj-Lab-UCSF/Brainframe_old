function brainframe(input_struct)

%Figure call
figure;

%Voxel or region resolution flag
voxUreg = input_struct.voxUreg;

%Connectome visualization binary flag
iscon = input_struct.iscon;

%Base brain atlas prep
brainat = input_struct.brain_atlas;
CCFinds = brainat;
CCFbin = brainat;
CCFbin(CCFbin > 0) = 1;
isomap = CCFbin;
isomap = smooth3(isomap,'box',5);
isobin = isomap;
isobin(isobin > 0) = 1;
CCFinds = CCFinds .* isobin;
surf1 = isosurface(isomap);
p1 = patch(surf1);
% v = get(p1);
% w = isonormals(isomap,p1);
if strcmp(input_struct.bgcolor,'w')
    set(p1,'FaceColor',[0.9 0.85 0.85],'EdgeColor','none','FaceAlpha',0.05);
elseif strcmp(input_struct.bgcolor,'k')% set the color, mesh and transparency level of the surface
    set(p1,'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.05);
else
    set(p1,'FaceColor',[0.95 0.9 0.9],'EdgeColor','none','FaceAlpha',0.0625);
end
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud
hold on

%Binary flag for per voxel (0) or per region (1) rendering
if voxUreg
    %Per-region rendering
    %Setting up data, per region
    groupid = input_struct.region_groups;
    idlist = unique(groupid);
    is0 = (idlist==0);
    idlist(is0) = [];
    scalefac = input_struct.xfac;
    cmap = input_struct.cmap;
    scalevec = input_struct.data;
    normvec = scalevec / mean(scalevec);
    normvec = normvec * scalefac;
    ptsz = input_struct.pointsize;
    sphpts = input_struct.sphere_npts;
    centroids = zeros(length(groupid),3);
    
    %Grouping regions for plotting
    centdex = 0;
    for i = 1:length(idlist)
        curid = idlist(i);
        reginds = find(groupid==curid);
        centers = [0 0 0];
        normvals = normvec(reginds);
        
        for k = 1:length(reginds)
            curreg = reginds(k);
            pcinds = find(CCFinds==curreg);
            [x,y,z] = ind2sub(size(CCFbin),pcinds);
            centroid = [mean(y) mean(x) mean(z)];
            if ~(isempty(input_struct.conmat)) && size(input_struct.conmat,1) == size(input_struct.data,1)
            	centroids(reginds(k),:) = centroid;
            end
            if scalevec(reginds(k)) > 0
                %Finding region centers, placing a sphere around them
                if input_struct.sphere
                    [sphx,sphz,sphy] = sphere(ceil(sphpts));
                    centroid = repmat(centroid,ceil(sphpts+1)^2,1) +...
                        (normvals(k))*[sphx(:) sphy(:) sphz(:)];
%                     [sphx,sphz,sphy] = sphere(ceil(normvals(k)));
%                     centroid = repmat(centroid,ceil(normvals(k)+1)^2,1) +...
%                         (normvals(k))*[sphx(:) sphy(:) sphz(:)];
                    centroid = repmat(centroid,sphpts,1);
                %Setting up diffuse random point clouds per region
                elseif ~input_struct.sphere && ~input_struct.centered(1)
                    rng(k);
                    randinds = randi(length(pcinds),ceil(normvals(k)*length(pcinds)^(1/3)),1);
%                     randinds = randi(length(pcinds),ceil(normvals(k)),1);
                    chosevox = [y(randinds) x(randinds) z(randinds)];
                    centroid = chosevox + rand(size(chosevox,1),3);
                %Creating random point clouds biased towards region centers
                elseif ~input_struct.sphere && input_struct.centered(1)
                    if length(input_struct.centered) > 1
                        centstrength = input_struct.centered(2);
                    else
                        centstrength = 0.5;
                    end
                    voxdists = ((y-centroid(1)).^2 + (x-centroid(2)).^2 + (z-centroid(3)).^2).^(centstrength);
                    dieinds = WeightedDie(voxdists);
%                     npts = ceil(normvals(k));
                    if mean(nonzeros(scalevec)) == 1
%                         npts = ceil(normvals(k)*length(pcinds));
                          npts = ceil(normvals(k));
                    else
                        npts = ceil(normvals(k)*length(pcinds)^(1/3));
                    end
                    rng(k);
                    choseinds = randi(length(dieinds),npts,1);
                    choseinds = dieinds(choseinds);
                    choseinds = pcinds(choseinds);
                    [x_chose,y_chose,z_chose] = ind2sub(size(CCFbin),choseinds);
                    chosevox = [y_chose x_chose z_chose];
                    centroid = chosevox + rand(size(chosevox,1),3); %.* multipliers;
                end
                centers = [centers;centroid];
            end
            centdex = centdex + 1;
        end
        
        %Plotting regions at center, rendered as sphere
        centers(1,:) = [];
        dists = sqrt(centers(:,1).^2 + centers(:,2).^2 + centers(:,3).^2);
        intensities = linspace(0.75,1,size(centers,1))';
        [~,sortinds] = sort(dists);          
        intensities = intensities(sortinds);
        ptcloud = pointCloud(centers,'Color',...
            repmat(cmap(i,:),size(centers,1),1),...
            'Intensity',intensities);
        pcshow(ptcloud,'MarkerSize',ptsz); hold on; 
    end
    
    %Rendering per-region connectivity map
    if iscon
        conmat = input_struct.conmat;
%         conmaxnorm = conmat / max(conmat(:));
        normmat = conmat / max(conmat(:));
        conscale = input_struct.con_rescale;
        conarch = input_struct.con_arch;
        if conarch > 1
            conarch = 1;
        elseif conarch < 0
            conarch = 0;
        end
        conwidth = input_struct.con_width;
        colormatrix = input_struct.con_cmap;
        congroups = input_struct.con_regiongroups;
        arrowWL = input_struct.conarrow_WL;
%         centers_ell = mean(centroids);
        conx_up = []; cony_up = []; conz_up = [];
        conx_dn = []; cony_dn = []; conz_dn = [];
        arrow_up1 = []; arrow_up2 = [];
        arrow_dn1 = []; arrow_dn2 = [];
        colors_up = []; colors_dn = [];
        arcol_up = []; arcol_dn = [];
        crossdir_up = []; crossdir_dn = [];
        tipang_up = []; tipang_dn = [];
%         u_up = []; v_up = []; w_up = [];
%         u_dn = []; v_dn = []; w_dn = [];
        for g = 1:size(conmat,1)
%             curcon = conmat(g,:);
%             curcon = curcon.';
%             curcon = floor(curcon*conscale);
%             normcon = normmat(g,:);
%             curcmap = concmap(g,:);
%             conx = []; cony = []; conz = [];
            for h = g:size(conmat,1)
                if (conmat(g,h) > 0 || conmat(h,g) > 0) && ~(h==g) && centroids(g,1) > 0 && centroids(h,1) > 0
                    curx = [centroids(g,2) centroids(h,2)].';
                    cury = [centroids(g,1) centroids(h,1)].';
                    curz = [centroids(g,3) centroids(h,3)].';

                    % calculate the direction of the connecting line between the two points
                    v = [curx(2)-curx(1); cury(2)-cury(1); curz(2)-curz(1)];
                    v_ = v/(v.' * v)^0.5;
                    
                    % calculate the rotation matrix to transform [1 0 0] to v
                    % (https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/897677#897677)
                    e1 = [1 0 0].';
                    G = [[e1.'*v_, -norm(cross(e1,v_),2), 0];...
                    [norm(cross(e1,v_),2), e1.'*v_, 0];...
                    [0, 0, 1]];
                    v_rej = v_ - (e1.'*v_*e1); v_rej = v_rej/(v_rej.'*v_rej)^0.5;
                    F = [e1, v_rej, cross(v_,e1)]\eye(3);
                    R = (F\eye(3))*G*F;
                    
                    % calculate the rotation of [0 0 1] using R to find the proper
                    % perpendicular direction
%                     if curx(1) > mean(centroids,2)
%                         e3 = [0 0 1].';
%                     else
%                         e3 = [0 0 -1].';
%                     end
                    e3 = [0 0 1].';
                    u = R*e3;
                    w = cross(v_,u);
                    
                    % calculate an ellipse centered at the midpoint of the line connecting the
                    % two points with major and minor axes aligned to v_ and u, respectively
                    midx = mean(curx);
                    midy = mean(cury);
                    midz = mean(curz);
                    midpt = [midx;midy;midz];
                    distxy = ((curx(2)-curx(1))^2 + (cury(2)-cury(1))^2 + (curz(2)-curz(1))^2)^0.5; % = 2*length major axis
                    fac = conarch; % perpendicular axis scale factor
                    nt = 50; % define the ellipse; 50 is arbitrary
                    narcs_gh = 2 * ceil(conmat(g,h)*conscale) + 2; % number of random arcs
                    narcs_hg = 2 * ceil(conmat(h,g)*conscale) + 2;
                    theta_gh = log(1 + normmat(g,h)) * pi/12;
                    theta_hg = log(1 + normmat(h,g)) * pi/12; 
                    rng(0);
                    rands_gh = rand(1,narcs_gh);
                    rands_hg = rand(1,narcs_hg);
                    randtheta_gh = theta_gh*(2*rands_gh - 1);
                    randtheta_hg = theta_hg*(2*rands_hg - 1);
                    gh_on = (conmat(g,h) > 0);
                    hg_on = (conmat(h,g) > 0);
                    uws_gh = repmat(u,1,narcs_gh).*cos(repmat(randtheta_gh,length(u),1)) + ...
                            repmat(w,1,narcs_gh).*sin(repmat(randtheta_gh,length(w),1));
                    uws_hg = repmat(u,1,narcs_hg).*cos(repmat(randtheta_hg,length(u),1)) + ...
                            repmat(w,1,narcs_hg).*sin(repmat(randtheta_hg,length(w),1));
                    if v(3) > 0
                        t_up = reshape(linspace(0,pi,nt),1,1,nt); % full ellipse is defined from 0 to 2*pi
                        t_dn = reshape(linspace(pi,2*pi,nt),1,1,nt);
                    else
                        t_up = -reshape(linspace(0,pi,nt),1,1,nt); % full ellipse is defined from 0 to 2*pi
                        t_dn = -reshape(linspace(pi,2*pi,nt),1,1,nt);
                    end
                    ellt_up = hg_on*(repmat(midpt,1,narcs_hg,nt) ...
                                + 0.5*distxy*repmat(v_,1,narcs_hg,nt).*cos(repmat(t_up,length(v_),narcs_hg,1)) ...
                                + 0.5*fac*distxy*repmat(uws_hg,1,1,nt).*sin(repmat(t_up,length(u),narcs_hg,1))); % parametric equation
                    ellt_dn = gh_on*(repmat(midpt,1,narcs_gh,nt) ...
                                + 0.5*distxy*repmat(v_,1,narcs_gh,nt).*cos(repmat(t_dn,length(v_),narcs_gh,1)) ...
                                + 0.5*fac*distxy*repmat(uws_gh,1,1,nt).*sin(repmat(t_dn,length(u),narcs_gh,1)));
                    if hg_on
                        ellx_up = squeeze(ellt_up(1,:,1:floor(1*nt))).';
                        elly_up = squeeze(ellt_up(2,:,1:floor(1*nt))).';
                        ellz_up = squeeze(ellt_up(3,:,1:floor(1*nt))).';
                        conx_up = [conx_up ellx_up];
                        cony_up = [cony_up elly_up];
                        conz_up = [conz_up ellz_up];
%                         endpt_up1 = squeeze(ellt_up(:,:,floor(0.45*nt))).';
                        endpt_up1 = mean(squeeze(ellt_up(:,:,floor(0.65*nt))).',1);
                        arrow_up1 = [arrow_up1; endpt_up1];
%                         endpt_up2 = squeeze(ellt_up(:,:,ceil(0.55*nt))).';
                        endpt_up2 = mean(squeeze(ellt_up(:,:,ceil(0.7*nt))).',1);
                        arrow_up2 = [arrow_up2; endpt_up2];
                        colors_up = [colors_up; repmat(colormatrix(congroups(h),:),size(ellx_up,2),1)];
                        arcol_up = [arcol_up; colormatrix(congroups(h),:)];
%                         arcol_up = [arcol_up; repmat(colormatrix(congroups(h),:),size(endpt_up1,1),1)];
%                         crossdir_up = [crossdir_up;cross(w,(endpt_up2 - endpt_up1)/norm((endpt_up2 - endpt_up1),2))];
%                         tipang_up = [tipang_up;(25*(conmat(g,h)/max(conmat(:)))+20)];
                    end
                    if gh_on
                        ellx_dn = squeeze(ellt_dn(1,:,1:floor(1*nt))).';
                        elly_dn = squeeze(ellt_dn(2,:,1:floor(1*nt))).';
                        ellz_dn = squeeze(ellt_dn(3,:,1:floor(1*nt))).';
                        conx_dn = [conx_dn ellx_dn];
                        cony_dn = [cony_dn elly_dn];
                        conz_dn = [conz_dn ellz_dn];
                        endpt_dn1 = mean(squeeze(ellt_dn(:,:,floor(0.75*nt))).',1);
%                         endpt_dn1 = squeeze(ellt_dn(:,:,floor(0.45*nt))).';
                        arrow_dn1 = [arrow_dn1; endpt_dn1];
%                         endpt_dn2 = squeeze(ellt_dn(:,:,ceil(0.55*nt))).';
                        endpt_dn2 = mean(squeeze(ellt_dn(:,:,ceil(0.8*nt))).',1);
                        arrow_dn2 = [arrow_dn2; endpt_dn2];
                        colors_dn = [colors_dn; repmat(colormatrix(congroups(g),:),size(ellx_dn,2),1)];
%                         arcol_dn = [arcol_dn; repmat(colormatrix(congroups(g),:),size(endpt_dn1,1),1)];
                        arcol_dn = [arcol_dn; colormatrix(congroups(g),:)];
%                         crossdir_dn = [crossdir_dn;cross(w,(endpt_dn2 - endpt_dn1)/norm((endpt_dn2 - endpt_dn1),2))];
%                         tipang_dn = [tipang_dn;(25*(conmat(h,g)/max(conmat(:)))+20)];
                    end
                    
%                     u_up = [u_up (ellx_up(nt,:)-ellx_up(nt-1,:))/distxy];
%                     v_up = [v_up (elly_up(nt,:)-elly_up(nt-1,:))/distxy];
%                     w_up = [w_up (ellz_up(nt,:)-ellz_up(nt-1,:))/distxy];
%                     u_dn = [u_dn (ellx_dn(nt,:)-ellx_dn(nt-1,:))/distxy];
%                     v_dn = [v_dn (elly_dn(nt,:)-elly_dn(nt-1,:))/distxy];
%                     w_dn = [w_dn (ellz_dn(nt,:)-ellz_dn(nt-1,:))/distxy];


%                     plot3([centroids(g,1) centroids(h,1)],[centroids(g,2) centroids(h,2)],[centroids(g,3) centroids(h,3)],'Color',curcmap,'LineWidth',curcon(h)*conwidth);
                end
                
            end

        end
        colors_up = num2cell(colors_up,2); colors_dn = num2cell(colors_dn,2);
        arcol_up = num2cell(arcol_up,2); arcol_dn = num2cell(arcol_dn,2);
%     crossdir_up = crossdir_up.'; crossdir_dn = crossdir_dn.';
        if ~isempty(conx_up) 
            pup = plot3(cony_up,conx_up,conz_up,'LineWidth',conwidth); hold on;
            set(pup,{'color'},colors_up);
%         aup = arrow(arrow_up1(:,[2,1,3]),arrow_up2(:,[2,1,3]),'TipAngle',25);
%         aup = arrow(arrow_up1(:,[2,1,3]),arrow_up2(:,[2,1,3]),'CrossDir',crossdir_up(:,[2 1 3]),'TipAngle',45,'Length',5);
            aup = arrow3(arrow_up1(:,[2,1,3]),(arrow_up2(:,[2,1,3])-arrow_up1(:,[2,1,3])),'0',arrowWL(1),arrowWL(2),'cone'); hold on;
            set(aup,{'FaceColor'},arcol_up); set(aup,{'EdgeColor'},arcol_up);
%         hq = quiver3(cony_up(end-1,:),conx_up(end-1,:),conz_up(end-1,:),v_up,u_up,w_up,...
%                 'LineWidth',1,'Color',[0 1 0]); hold on;
% %         hq.NodeChildren(2).Visible = 'off';
        end
        if ~isempty(conx_dn)
            pdn = plot3(cony_dn,conx_dn,conz_dn,'LineWidth',conwidth); hold on;
            set(pdn,{'color'},colors_dn);
%         adn = arrow(arrow_dn1(:,[2,1,3]),arrow_dn2(:,[2,1,3]),'TipAngle',25);
%         adn = arrow(arrow_dn1(:,[2,1,3]),arrow_dn2(:,[2,1,3]),'CrossDir',crossdir_dn(:,[2 1 3]),'TipAngle',45,'Length',5);
            adn = arrow3(arrow_dn1(:,[2,1,3]),(arrow_dn2(:,[2,1,3])-arrow_dn1(:,[2,1,3])),'0',arrowWL(1),arrowWL(2),'cone'); hold on;
            set(adn,{'FaceColor'},arcol_dn); set(adn,{'EdgeColor'},arcol_dn);
%         hq = quiver3(cony_dn(end-1,:),conx_dn(end-1,:),conz_dn(end-1,:),v_dn,u_dn,w_dn,...
%                 'LineWidth',1,'Color',[1 0 0]); hold on;
%         hq.NodeChildren(2).Visible = 'off';
% %                 plot3(cony.',conx.',conz.','Color',curcmap,'LineWidth',conwidth); hold on;
        end
    end
else
    %Per-voxel rendering
    %Data mapping and prep, per voxel
    interpmap = input_struct.data;
    if size(interpmap,2) <= 1
        interpmap = reshape(interpmap,size(CCFbin));
    end    
    testvals = interpmap(:);
    sumvals = sum(testvals);
    sortvals = sort(testvals,'descend');
    cumsumvals = cumsum(sortvals);
    fracvals = cumsumvals / sumvals;
    interpmap = CCFbin .* interpmap;
    cmap = input_struct.cmap;
    interpmap_f = smooth3(interpmap,'box',[7 1 1]);
    interpmap_f = CCFbin .* interpmap_f;
    
    %Setting up binning and misc for visualization
    nbin = input_struct.nbin;
    rangei = linspace(0,input_struct.voxthresh,nbin+1);
    rangei = fliplr(rangei);
    xfac = input_struct.xfac;
    ptsz = input_struct.pointsize;
    mins = sortvals(fracvals<rangei(1) & fracvals>=rangei(1+1));
    meanmins = mean(mins);
    
    %Looping through bins and visualizing point clouds
    for i = 1:nbin
        bounds = sortvals(fracvals<rangei(i) & fracvals>=rangei(i+1));
        bound1 = bounds(1); bound2 = bounds(end);
        [x,y,z] = ind2sub(size(interpmap_f),find(interpmap_f>=bound2 & interpmap_f<bound1));
        xyzmap = [y x z];
        ifac = xfac * (mean(bounds)/meanmins);
        xyz_jitter = repmat(xyzmap,floor(ifac)+1,1) + rand(size(xyzmap,1)*(1+floor(ifac)),size(xyzmap,2));
        ptcloud = pointCloud(xyz_jitter,'Color',repmat(cmap(i,:),size(xyz_jitter,1),1),...
            'Intensity',repmat(0.1*i,size(xyz_jitter,1),1));
        pcshow(ptcloud,'MarkerSize',ptsz); hold on;
        clear x y z xyzmap ptcloud xyz_jitter
    end
end    

%Changing axis properties for visualization
set(gcf,'color',input_struct.bgcolor);
set(gca,'color',input_struct.bgcolor);
ax = gca;
set(ax,'XColor','none','YColor','none','ZColor','none');

%Saving and closing or opening .fig GUI
savenclose = input_struct.savenclose;
imglab = input_struct.img_labels;
imgtype = input_struct.img_format;
set(gcf, 'InvertHardCopy', 'off'); 
if savenclose
    view([0, 0, 1]);
    ax = gca;
    set(ax,'XColor','none','YColor','none','ZColor','none');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    saveas(gcf,[imglab '_sagittal'],imgtype);
    
    view([-1, 0, 0]);
    ax = gca;
    set(ax,'XColor','none','YColor','none','ZColor','none');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    saveas(gcf,[imglab '_axial'],imgtype);
    
    view([0, -1, 0]);
    ax = gca;
    set(ax,'XColor','none','YColor','none','ZColor','none');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    saveas(gcf,[imglab '_coronal'],imgtype);
    close
    clear isomap
end

end

function [weighted_die] = WeightedDie(distarr)
% This function creates a weighted die of indices based on an input of
% distances, where the "weight" is inversely proportional to the distance.

% 1. Create a "similarity" capped at 100 with a minimum of 1
distarr = reshape(distarr,[length(distarr),1]);
sim = ones(length(distarr),1) ./ distarr;
sim_n = sim / min(sim);
sim_n = (5000/max(sim_n))*sim_n;
sim_r = ceil(sim_n);

weighted_die = [];
for i = 1:length(sim_r)
    weighted_die = [weighted_die i*ones(1,sim_r(i))];
end
weighted_die = weighted_die.';

end









