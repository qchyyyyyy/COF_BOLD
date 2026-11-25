function [z_stats,p_values] = fit_vertex_glm(X,Y)
lm = fitlm(X, Y);
b = lm.Coefficients.Estimate;
se_b2 = lm.Coefficients.SE(2);
z_stats = b(2) / se_b2;
p_values = 2 * (1 - tcdf(abs(b(2)) / se_b2, lm.DFE));