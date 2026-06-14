%lambda:nutrient penetration length,
%P_res:mechanical resistance exerted by neighbouring cells
%Therapies: 2-element vector, [K_RT,K_CT]
function metrics=CA_Avascular_Tumors(lambda,P_res,Therapies)

radiotherapy=0;
chemotherapy=0;

if nargin==3
    RT=Therapies(1);
    QT=Therapies(2);
    if RT~=0
        radiotherapy=1;
    end
    if QT~=0
        chemotherapy=1;
        Dl=0.8*lambda;
    end

end

margin=5;

%Treatment

t_RT=50;
t_QT=50;


%Cell states

empty=uint8(0);
proliferating=uint8(1);
quiescent=uint8(2);
nectrotic=uint8(3);


%Parameters

N=320;  %Grid length
n0=1; %Nutrients available at the EM
P_div=0.5; %Probability of Division
P_death=0.005; %Probability of spontaneous death (apoctosis?)
%P_res=0.2; %Pressure resistance exerted by the neighbouring cells
%lambda=10; %Diffusion characteristic length sqrt(D/k)
T_max=12; %Number of iterations in necrotic state before the cell dies


%Input Parameters
%{
lambda = input('Lambda: 3');
P_res = input('P_res: 0.95');

%}

%Neighborhood filter
Neighbors8 = ones(3,3,'double');
Neighbors8(2,2) = 0;


finish=0;

%% Initialization

T=ceil((N/2)/(P_div))+100+100*(radiotherapy+chemotherapy);
%time=1:T;
grid_states=zeros(N,N,'uint8');
grid_states(round(N/2),round(N/2))=proliferating;
new_cells=zeros(N,N,'logical');

tumor_mask=logical(grid_states~=empty);

quiescency_counter=zeros(N,N,'uint8');


Fd_interface = nan(1,T,'single');
Fd_living_cells = nan(1,T,'single');
gamma = nan(1,T,'single');
gamma2 = nan(1,T,'single');

S = nan(1,T,'single');
V = nan(1,T,'single');
V2 = nan(1,T,'single');

Nf=nan(1,T,'single');
Lf=nan(1,T,'single');
Pf=nan(1,T,'single');

i=0;
%%
max_reps=10;
for m=1:max_reps
    i=0;
while finish==0 && any(grid_states(:)==proliferating|grid_states(:)==quiescent)
    i=i+1;

    %Add new cells to the tumor
    tumor_mask=tumor_mask|new_cells;
    new_cells=zeros(N,N,'logical');

    %Compute neighbours density field
    rho=imfilter((single(tumor_mask)),Neighbors8,0);

    %Compute nutrient field
    grid_nutrients=get_nutrient_field(tumor_mask,lambda);

    %Cells die with P_death
    grid_states(P_death>rand(N,N)&tumor_mask>0)=nectrotic;

    %Living cells can proliferate or be quiescent
    living_cells=(tumor_mask>0&grid_states~=nectrotic);

    proliferation=grid_nutrients>rand(N,N)&living_cells;
    quiescency=proliferation==0&living_cells;
    
    grid_states(proliferation)=proliferating;
    grid_states(quiescency)=quiescent;

    %Update quiescency counter & cells with count>Tmax become necrotic
    quiescency_counter(proliferation)=0;
    quiescency_counter(quiescency)=quiescency_counter(quiescency)+1;
    grid_states(quiescency_counter>T_max)=nectrotic;

    %Apply therapy
    if radiotherapy==1 && i>t_RT
        grid_states(P_RT(grid_nutrients,grid_states,RT)>rand(N,N)&grid_states~=0)=nectrotic;
    end
    if chemotherapy==1 && i>t_QT
        grid_states(P_QT(grid_states,QT,Dl)>rand(N,N)&grid_states==1)=nectrotic;
    end



    %Compute the probability of proliferating
    P = P_div*(max(0,1-P_res.*rho/8));%.*grid_nutrients;
    P(grid_states~=proliferating)=0;
    duplicate=find(P>rand(N,N));

    %Loop over the proliferating cells in a random order
    duplicate=duplicate(randperm(numel(duplicate)));
    all_directions=[-1-1i,-1,-1+1i,-1i,1i,1-1i,1,1+1i];
    for d=duplicate'
        [x,y]=ind2sub([N,N], d);
        directions=[];
        sites=0;
        q=0;
        for l=-1:1
            for j=-1:1
                if (l ~= 0 || j ~= 0) 
                    q=q+1;
                    if ((x+l)<=N && (y+j)<=N) && ((x+l)>0 && (y+j)>0) && tumor_mask(x+l,y+j)==0 && new_cells(x+l,y+j)==0 
                        sites=sites+1;
                        directions=[directions,all_directions(q)];
                    end
                end
            end
        end
        if sites>0
            direction=directions(randi(sites));
            new_cells(x+real(direction),y+imag(direction))=new_cells(x+real(direction),y+imag(direction))+1;
        end
        
    end

    if all([tumor_mask(margin,:),tumor_mask(end-margin,:),tumor_mask(:,margin)',tumor_mask(:,end-margin)']==0)
        finish=0;
    else
        finish=1;
    end


    %Compute metrics
        %Fractal dimension
    empty_neighbors = imfilter(single(~tumor_mask), Neighbors8, 0) > 0;
    interface = living_cells & empty_neighbors;

    Fd_interface(i)=box_counting(interface);
    Fd_living_cells(i)=box_counting(living_cells);

        %Order parameter (S~V^gamma)
    S(i)=nnz(interface);
    V(i)=nnz(living_cells);
    V2(i)=nnz(tumor_mask);
    R_exp=sqrt(V2(i)/pi);
    num_samp=100;
    % {

    if i>300
        num_samp=round(i/2);
    end
    %}
    if i>num_samp
        gamma(i)=get_gamma(S(i-num_samp:i),V(i-num_samp:i));
    end
    if i>num_samp
        gamma2(i)=get_gamma(S(i-num_samp:i),V2(i-num_samp:i));
    end

        %Real vs expected state cells ratios
    Nf(i)=get_necrotic_ratio(grid_states,R_exp,lambda,nectrotic,T_max);
    Lf(i)=get_living_ratio(living_cells,R_exp,lambda,T_max);
    Pf(i)=get_proliferating_ratio(grid_states,R_exp,lambda,proliferating,T_max);
    
    q=i;
end

if q >= 20
    break
else
    grid_states(round(N/2),round(N/2))=proliferating;
    new_cells=zeros(N,N,'logical');
    tumor_mask=logical(grid_states~=empty);
    quiescency_counter=zeros(N,N,'uint8');

    if m==max_reps
        metrics.fail=1;
    end
end

end

hold off

value1=int64(max(q*4/5,1));
value2=q;


metrics.Fd_living_cells=mean(Fd_living_cells(value1:value2), 'omitnan');
metrics.Fd_interface=mean(Fd_interface(value1:value2), 'omitnan');
metrics.gamma=mean(gamma(value1:value2), 'omitnan');
metrics.gamma2=mean(gamma2(value1:value2), 'omitnan');
metrics.Nf=mean(Nf(value1:value2), 'omitnan');
metrics.Lf=mean(Lf(value1:value2), 'omitnan');
metrics.Pf=mean(Pf(value1:value2), 'omitnan');

if finish==0
    metrics.LivingMassRatio=nnz(grid_states(grid_states==proliferating|grid_states==quiescent))/nnz(grid_states);
    metrics.death_iterations=q;
    metrics.Mass=nnz(grid_states);
else
    metrics.LivingMassRatio=nnz(grid_states(grid_states==proliferating|grid_states==quiescent))/nnz(grid_states);
    metrics.Mass=nnz(grid_states);
    
end


% Plot cell states only (opcional)
% {
plot_states(grid_states)
%}

end


%
%Functions definition

    %Compute slope for gamma

function gamma_fit = get_gamma(S_samples, V_samples)

    S_samples = S_samples(:);
    V_samples = V_samples(:);

    valid = isfinite(S_samples) & isfinite(V_samples) & ...
            S_samples > 1 & V_samples > 1;

    S_samples = S_samples(valid);
    V_samples = V_samples(valid);

    if numel(S_samples) < 3
        gamma_fit = NaN;
        return
    end

    x = log(V_samples);
    y = log(S_samples);

    p = polyfit(x, y, 1);

    gamma_fit = p(1);

end

%Therapies probabilities

function p_rt=P_RT(nutrient_grid, grid, RT)

    p_rt=nutrient_grid*RT;
    p_rt(grid==0)=0;
end

function p_qt=P_QT(grid,QT,Dl)
    d=bwdist(~grid);
    C=exp(-d./Dl);
    p_qt=QT*C;
    p_qt(grid==0)=0;

end

%Nutrient map

function grid_nutrients=get_nutrient_field(tumor_mask,lambda)

    dist=single(bwdist(~tumor_mask));
    grid_nutrients=exp(-dist./lambda);
end

%fractal dimension

function Fd=box_counting(mask)

    mask = logical(mask);

    if nnz(mask) < 20
        Fd = NaN;
        return
    end

    [rows, cols] = find(mask);

    if isempty(rows)
        Fd = NaN;
        return
    end

    sizes = [2 4 8 16 32];

    if numel(sizes) < 3
        Fd = NaN;
        return
    end

    counts = zeros(size(sizes));

    [N,M] = size(mask);

    for k = 1:numel(sizes)
        s = sizes(k);

        dimN = floor(N/s) * s;
        dimM = floor(M/s) * s;

        cropped = mask(1:dimN, 1:dimM);

        blocks = reshape(cropped, s, dimN/s, s, dimM/s);
        blocks = permute(blocks, [2 4 1 3]);

        occupied_blocks = any(blocks, [3 4]);
        counts(k) = sum(occupied_blocks(:));
    end

    valid = counts > 0;

    if nnz(valid) < 3
        Fd = NaN;
        return
    end

    p = polyfit(log(1 ./ sizes(valid)), log(counts(valid)), 1);
    Fd = p(1);

end

%real/expected cells:

function  Nf=get_necrotic_ratio(states,R,lambda,necrotic,T_max)
    
    necrotic_cells=sum(states(:)==necrotic);
    ringN=(log(T_max)-1/lambda)*lambda;

    if R<20 || R<ringN*1.5
        Nf=NaN;
        return
    end

    expected_necrotic_cells=pi*(R-ringN)^2;
    Nf=necrotic_cells/expected_necrotic_cells;
end


function  Pf=get_proliferating_ratio(states,R,lambda,proliferating,Tmax)
    if R<20
        Pf=NaN;
        return
    end
    proliferating_cells=sum(states(:)==proliferating);

    xf=log(Tmax)*lambda;
    expected_proliferating_cells = 2 * pi * R * lambda * (exp(-1/lambda) - exp(-xf/lambda));

    Pf=proliferating_cells/expected_proliferating_cells;
end

function  Lf=get_living_ratio(living_mask,R,lambda,T_max)
    if R<20
        Lf=NaN;
        return
    end
    living_cells=sum(living_mask(:));
    ringL=(log(T_max)-1/lambda)*lambda;
    expected_living_cells=2*pi*(R-ringL/2)*ringL;%*(1-log(2));
    Lf=living_cells/expected_living_cells;
end

%final plot of cell states

function plot_states(grid_states)


    figure('Position', [100, 100, 600, 600], 'Color', 'w')

    colors=[1,1,1;0,0.8,0;1,0,0;0,0,0];
    
    ax = axes;
    imshow(grid_states, [0 3]);
    colormap(ax, colors);
    axis(ax, 'equal', 'tight');
    hold(ax, 'on');
    

    ax.Position = [0.12, 0.16, 0.76, 0.68];
    
    h0 = plot(ax, nan, nan, 's', ...
        'MarkerFaceColor', colors(1,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 15);
    
    h1 = plot(ax, nan, nan, 's', ...
        'MarkerFaceColor', colors(2,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 15);
    
    h2 = plot(ax, nan, nan, 's', ...
        'MarkerFaceColor', colors(3,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 15);
    
    h3 = plot(ax, nan, nan, 's', ...
        'MarkerFaceColor', colors(4,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 15);
    
    lgd = legend(ax, [h0, h1, h2, h3], ...
        {'ECM', 'Proliferating', 'Quiescent', 'Necrotic'}, ...
        'Orientation', 'horizontal', ...
        'Location', 'southoutside');
    
    legend(ax, 'boxoff');
    lgd.FontSize = 19;
    
    drawnow;
    
    center_plot = ax.Position(1) + ax.Position(3)/2;
    lgd.Position(1) = center_plot - lgd.Position(3)/2;
    lgd.Position(2) = 0.03;
    
    hold(ax, 'off');
end
