%edited on 02-27-14 by Shashaank Vattikuti

function f=calc_freq(S,dt)
%(S)pikes is assumed to be a 2d grid (neuron index (rows) x iteration (columns)), where
%1=spike and 0=no spike
%frequency is the frequency of spikes in a 100 msec window; samples are
%non-overlapping
%code could be modified for spatial and temporal overlap

    iterations=size(S,2);
    win_sec=0.1; %sampling window in seconds
    win=win_sec*1000/dt; %sampling window in iterations
    k=0;
    for i=1:win:iterations
        k=k+1;
        s=sum(S(:,(k-1)*win+1:k*win),2);
        f(:,k)=s/win_sec; %calculate Hz
    end
                
end