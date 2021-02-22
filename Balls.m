%% Two-dimensional unit balls

[x,y] = meshgrid(linspace(-1,1,200));

large_p = 2.^(1:10);

for k = 1: 10
    p = large_p(k);
    contour(x,y,abs(x).^p + abs(y).^p,'LevelList',1, 'linewidth', 3);
    title(['p = ', num2str(large_p(k))],'FontSize',14);
    axis equal
    axis([-1.5 1.5 -1.5 1.5])
    pause
end

%%
small_p = 2:-0.1:1;

for k = 1:11
    p = small_p(k);
    contour(x,y,abs(x).^p + abs(y).^p,'LevelList',1, 'linewidth', 3);
    title(['p = ', num2str(small_p(k))],'FontSize',14);
    axis equal
    axis([-1.5 1.5 -1.5 1.5])
    pause
end