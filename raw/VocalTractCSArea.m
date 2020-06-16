function ca = VocalTractCSArea(posterior, anterior,dims,points, make_plot)
% VOCALTRACTCSAREA - Derives cross-sectional area function for vocal tract
%
% SYNOPSIS: ca = VocalTractCSArea(posterior, anterior, size, points)
%
% INPUT: posterior: n x 2 integer array of x,y points corresopnding to a
% trace ofthe vocal tract from the glottis to the lips along the back
% pharyngeal wall
%        anterior: n x 2 integer array of x,y points corresopnding to a
% trace of the vocal tract from the glottis to the lips along the front
% pharyngeal wall
%        dims: integer corresponding to the dimensions of the MRI frame
%        points: number of evenly spaced points along the vocal tract for
%        which cross-sectional area should be calculated
%        make_plot: boolean indicating whether to show series of plots
%        visualizing the function results
%
% OUTPUT: ca: vector of length 'points', where each element contains the
% pixelwise area corresponding to each point along the vocal tract
%
% TO-DO: -Function will act strangely if glottis is vertical (or if
% mouth is horizontal)
%        -Runs in ~1.94s/frame; ways to speed up?
np = size(posterior,1);

%Constrain boundaries to integer values
posterior = round(posterior);
anterior = round(anterior);



%%Make the horizontal end of the posterior trace closer to 1
posdif = diff(posterior);
if(mean(abs(posdif(round(np/2):np-1,1))) > mean(abs(posdif(round(np/2):np-1,2))))
    posterior = flipud(posterior);
end

%%Order points so anterior and posterior correspond
normal= norm(posterior(1,:)-anterior(1,:));
opposite= norm(posterior(1,:) - anterior(size(anterior,1),:));

if opposite < normal
    anterior = flipud(anterior);
end

%%Ensure that extra mouth projection is perfectly horizontal

%Define proper horizontal direction
if(posterior(1,1) < mean(posterior(:,1)))
    extrx = min(anterior(1,1),posterior(1,1));
else 
    extrx = max(anterior(1,1),posterior(1,1));
end
   
%Define proper vertical direction
if(posterior(np,2) < mean(posterior(:,2)))
    extry = min(anterior(np,2),posterior(np,2));
else 
    extry = max(anterior(np,2),posterior(np,2));
end

posterior(1,1) = extrx;
anterior(1,1) = extrx;
anterior(np,2) = extry;
posterior(np,2) = extry;

%posterior
%anterior
%    plot(posterior(:,1),posterior(:,2),anterior(:,1),anterior(:,2));
%    pause();

% Info related to mouth opening
mouthm = (posterior(1,2)-anterior(1,2))/(posterior(1,1)-anterior(1,1)); %Slope of Line Passing Through Mouth

%%If Lips are Vertical
if abs(mouthm) == Inf
    mouthb = posterior(1,1);
else
    mouthb = posterior(1,2) - mouthm*posterior(1,1);
end

%%Intercepts of Upper and Lower Lines Perpendicular to Mouth
upperb = posterior(1,2) + (1/mouthm)*posterior(1,1);
lowerb = anterior(1,2) + (1/mouthm)*anterior(1,1);

airway = airwayPrep(posterior, anterior, dims, [mouthm, mouthb, upperb, lowerb]);
extra = double(airway~=-4);
contoured = regionGrowth({posterior(np,:), anterior(np,:)},...
    10,airway);
cent_curve = calculateCentroid(contoured, [mouthm, mouthb, upperb, lowerb]);
disp(cent_curve)
ca = findArea(cent_curve, posterior, anterior, points);

ca(:,8) = (pi*(ca(:,8)/2).^2)/(290); %2.9mm2/pixel in plane resolution (Sorenson et al 2017) * mm2 to cm2 conversion

%Plotting%
if(make_plot)
    subplot(2,2,1)
    plot(posterior(:,1),posterior(:,2),anterior(:,1),anterior(:,2));
    title('Air/Tissue Boundary of Vocal Tract');
    xlabel('X-Wise Pixel Value')
    ylabel('Y-Wise Pixel Value')
    set(gca, 'YDir','reverse');
    subplot(2,2,2)
    
    heatmap((contoured.*extra)', 'GridVisible', 'off','Colormap', parula);
    title('Airway Heatmap (Wamer Colors are Further from the Glottis)');
    xlabel('X-Wise Pixel Value')
    ylabel('Y-Wise Pixel Value')
    
    subplot(2,2,3)
    plot(posterior(:,1),posterior(:,2),anterior(:,1),anterior(:,2), cent_curve(:,2),cent_curve(:,3));
    hold on
    for i = 1:size(ca,1);
        plot(ca(i,[4 6]), ca(i,[5 7]),'-');
    end
    hold off
    title('Distances Calculated Along the Vocal Tract')
    xlabel('X-Wise Pixel Value')
    ylabel('Y-Wise Pixel Value')
    set(gca, 'YDir','reverse');
    subplot(2,2,4)
    plot(ca(:,1),ca(:,8));
    title('Cross Sectional Area Function')
    xlabel('Distance from Lips')
    ylabel('Cross-Sectional Area (cm^2)')
end



function ap = airwayPrep(pcontour, acontour, dims, mouth_info)
mouthm = mouth_info(1);
mouthb = mouth_info(2);
upperb = mouth_info(3);
lowerb = mouth_info(4);

disp(dims)
[x y] = find(zeros(dims(1),dims(2))+1);
poly = [pcontour ; flipud(acontour)];
px = poly(:,1);
py = poly(:,2);
[in on] = inpolygon(x,y,px,py);
ap = zeros(dims(1),dims(2)) -3;

for i=1:size(x)
    if abs(mouthm) ~= Inf
        btest = y(i) + (1/mouthm)*x(i);
        if btest >= upperb && btest <= lowerb && y(i) - mouthm*x(i) - mouthb > 0
            ap(x(i),y(i)) = -4;
        end
    else
        if y(i) >= upperb && y(i) <= lowerb && x(i) < mouthb
            ap(x(i),y(i)) = -4;
        end
    end
end

ap(sub2ind(size(ap),x(in),y(in))) = -2;


function region = regionGrowth(seed,seed_value,map)
% Non-negative integers: Airway space that has been traversed
% -1: Airway space in queue for traversing
% -2: Airway space that has never been encountered
% -3: Space ineligible for traversing

extra = double(map ~= -4);

remaining = seed;
temp = {};
region = map;

for a=1:length(seed)
    spoint = seed{a};
    region(spoint(1), spoint(2)) = seed_value;
end

first=seed{1};
last= seed{length(seed)};

m = (last(2)-first(2))/(last(1) - first(1));
b = first(2) - m*first(1);


col = repmat(1:20,20,1);
row = col';

while length(remaining)>0
    %Pop next element from cell array
    point = remaining{1};
    remaining(1) = [];
    
    %Check to see which neighboring points exist.
    if point(1) < size(region,1)
        temp{length(temp)+1}=[point(1)+1, point(2)];
    end
    if point(1) > 1
        temp{length(temp)+1}=[point(1)-1, point(2)];
    end
    if point(2) < size(region,2)
        temp{length(temp)+1}=[point(1), point(2)+1];
    end
    if point(2) > 1
        temp{length(temp)+1}=[point(1), point(2)-1];
    end
    
    %Of the neighboring points, determine which neighbor determines the point's value
    if(region(point(1), point(2))) == -1
        windex = min(find(cellfun(@(x) region(x(1),x(2)) > -1, temp)));
        if length(windex) > 0
            winner = temp{windex};
            if b + m*point(1) - point(2) > 0 %This may cause an issue if the VT is not oriented correctly
                region(point(1), point(2)) = region(winner(1), winner(2)) + 1;
            elseif b + m*point(1) - point(2) == 0
                region(point(1), point(2)) = region(winner(1), winner(2)) - 1;
            else
                region(point(1), point(2)) = region(winner(1), winner(2)) -1;
            end
        end
    end
    
    %And which points still need to be traversed
    newpoints = temp(find(cellfun(@(x) region(x(1), x(2)) == -2 | region(x(1), x(2)) == -4, temp)));
    remaining = horzcat(remaining, newpoints);
    %Mark that it has been traversed
    for a=1:length(newpoints)
        npoint = newpoints{a};
        region(npoint(1), npoint(2)) = -1;
    end
    
    temp = {};
end


function centcurve = calculateCentroid(array, mouth_info)
mouthm = mouth_info(1);
mouthb = mouth_info(2);
upperb = mouth_info(3);
lowerb = mouth_info(4);


uni = unique(array);
suitable = uni(find(uni>=0))';

centcurve=[];


for a=suitable
    [row, col] = find(array==a);
    avgx = mean(row);
    avgy = mean(col);
    if abs(mouthm) ~= Inf
        btest = avgy + (1/mouthm)*avgx;
        if ~(btest >= upperb && btest <= lowerb && avgy - mouthm*avgx - mouthb > 0)
            centcurve = [centcurve;[a avgx avgy]];
        end
    else
        if ~(avgy >= upperb && avgy <= lowerb && avgx < mouthb)
            centcurve = [centcurve;[a avgx avgy]];
        end
    end
end


function fa = findArea(curve, posterior, anterior, points)
nearest = 20;
%Smoothed Curves
xsmooth = fit(curve(:,1), curve(:,2),'poly2');
ysmooth = fit(curve(:,1), curve(:,3),'poly2');


fa = [];
for i = 1:points
    a = max(curve(:,1))-(i-1)*(max(curve(:,1))-min(curve(:,1)))/(points+1);
    %Get line coefficients for point
    m=(ysmooth(a+1e-3)-ysmooth(a))/(xsmooth(a+1e-3)-xsmooth(a));
    b=ysmooth(a) - m*xsmooth(a);
    pb = ysmooth(a) + (1/m)*xsmooth(a);
    
    ant = knnsearch(anterior, [xsmooth(a) ysmooth(a)], 'K',nearest);
    pos = knnsearch(posterior, [xsmooth(a) ysmooth(a)], 'K', nearest);
    
    [dummy,ai] = min(abs(anterior(ant,2)+(1/m)*anterior(ant,1) - pb));
    [dummy,pi] = min(abs(posterior(pos,2)+(1/m)*posterior(pos,1) - pb));
    
    chosen = [anterior(ant(ai),:);posterior(pos(pi),:)];
    fa = [fa;[i, xsmooth(a), ysmooth(a), anterior(ant(ai),:), posterior(pos(pi),:),pdist(chosen)]];
    
end
