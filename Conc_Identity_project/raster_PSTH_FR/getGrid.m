function[gridX, gridY] = getGrid(numStim)
% Gives gridX and gridY if 20 stimuli (concentration) or 16 stimuli
% (identity)
    if numStim == 20
        gridX = 5;
        gridY = 4;
    else
        gridX = 4;
        gridY = 4;
    end
end