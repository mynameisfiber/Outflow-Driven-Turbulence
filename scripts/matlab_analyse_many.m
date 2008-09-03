size=120
buffer=3
hf= figure('visible','off');
for frame=250:250:5750
	file = ['output-' num2str(frame, '%.8d')]
	u = load(file);
   view = reshape(u(:,1),size+2*buffer,size+2*buffer,size+2*buffer);
	view = view(buffer+1:size+buffer,buffer+1:size+buffer,buffer+1:size+buffer);

	[status, result] = system(['grep ' num2str(frame, '%.10d') ' outputtimes | cut -d " " -f 4']);
	
	clf();
	h = pcolor(log10(abs(reshape(mean(view,2),size,size))));colorbar;
	set(h,'edgecolor','none');
	title(['<rho> at t=' result]);
	xlabel('x'); ylabel('z');
	saveas(hf,['images/pcolor_mean2_' num2str(dim) '_' num2str(frame, '%.5d') '.png']);

	clf();
	h = pcolor(log10(abs(reshape(mean(view,3),size,size))));colorbar;
	set(h,'edgecolor','none');
	title(['<rho> at t=' result]);
	xlabel('x'); ylabel('y');
	saveas(hf,['images/pcolor_mean3_' num2str(dim) '_' num2str(frame, '%.5d') '.png']);
	
	clf();
	h = pcolor(log10(abs(view(:,:,size/2))));colorbar;
	set(h,'edgecolor','none');
	title(['Rho slice at z=33 & t=' result]);
	xlabel('x'); ylabel('z');
	saveas(hf,['images/pcolor_slicemid_' num2str(dim) '_' num2str(frame, '%.5d') '.png']);
	
	clf();
	[x,y,z] = meshgrid(1:1:size,1:1:size,1:1:size);        
	%xslice = [size/2,size]; yslice = [size/2, size]; zslice = [1,size/2];
	xslice = [size]; yslice = [size]; zslice = [1];
	h = slice(x,y,z,view,xslice,yslice,zslice);
	set(h,'FaceColor','interp','EdgeColor','none')
	axis tight
	box on
	set(gca, 'CLim', [min(view(:)), max(view(:))]);
	colormap (jet(24))
	colorbar('horiz')
	set(gcf,'Renderer','OpenGL')
	%set(gcf,'Renderer','zbuffer')
	data = smooth3(view,'box',5);
	isoval = 1.001e-04;
	%isoval = 1.001e-4
	h = patch(isosurface(data),...
		'FaceColor','red',...
		'EdgeColor','none',...
		'AmbientStrength',.2,...
		'SpecularStrength',.7,...
		'DiffuseStrength',.4,...
		'FaceAlpha',.5);
	isonormals(data,h)
	%patch(isocaps(data),...
	%	'FaceColor','interp',...
	%	'EdgeColor','none')
	daspect([1,1,1])
	axis tight
	xlabel 'x'
	ylabel 'y'
	zlabel 'z'
	title(['Rho at t=' result])
	view(2);
	camlight right
	camlight left
	lighting phong
	saveas(hf,['images/isosurface_' num2str(frame, '%.8d') '.png'])
	
	clf();
	[x,y,z] = meshgrid(1:1:size,1:1:size,1:1:size);        
	%xslice = [size/2,size]; yslice = [size/2, size]; zslice = [1,size/2];
	xslice = [size]; yslice = [size]; zslice = [1];
	h = slice(x,y,z,view,xslice,yslice,zslice);
	set(h,'FaceColor','interp','EdgeColor','none')
	axis tight
	box on
	set(gca, 'CLim', [min(view(:)), max(view(:))]);
	colormap (jet(24))
	colorbar('horiz')
	set(gcf,'Renderer','OpenGL')
	%set(gcf,'Renderer','zbuffer')
	data = smooth3(view,'box',5);
	isoval = 1.001e-04;
	%isoval = 1.001e-4
	h = patch(isosurface(data,15),...
		'FaceColor','red',...
		'EdgeColor','none',...
		'AmbientStrength',.2,...
		'SpecularStrength',.7,...
		'DiffuseStrength',.4,...
		'FaceAlpha',.5);
	isonormals(data,h)
	%patch(isocaps(data),...
	%	'FaceColor','interp',...
	%	'EdgeColor','none')
	daspect([1,1,1])
	axis tight
	xlabel 'x'
	ylabel 'y'
	zlabel 'z'
	title(['Rho w/ isosurface=15 at t=' result])
	view(2);
	camlight right
	camlight left
	lighting phong
	saveas(hf,['images/isosurface_15_' num2str(frame, '%.8d') '.png'])
	
end

