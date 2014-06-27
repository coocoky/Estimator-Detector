function reportChange(sampleIndex,...
                      windowNumber,...
                      distance)

fprintf('Changed detected around sample %d\n', sampleIndex);
fprintf('Window pair that detected the change: window%d\n', windowNumber);
fprintf('Distance = %.6f\n', distance);
fprintf('\n');
end % function reportChange