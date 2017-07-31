function plotHandle = PlotHistory(x,t,U,p,parentHandle,idx)

   %% Position and eventually grab figure
   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end
   figure(parentHandle);

   %% Extract number of components
   numComp = size(idx,2);

   %% Grid
   [T,X] = meshgrid(t,x);

   %% Plots
   for k = 1:numComp
     subplot(1,numComp,k)
     % pcolor(X,T,U(:,idx(:,k))'); shading interp; view([0 90]);
     surf(X,T,U(:,idx(:,k))'); shading interp; view([0 90]);
     title(['U_' num2str(k)]);
   end

   %% Save
   print -dtiff history.tiff

end

