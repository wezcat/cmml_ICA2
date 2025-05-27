# 组织成数据框
df <- data.frame(
  row.names = c("origin_data", "seurat", "scanvi", "scpoli"),
  ASW_label = c(0.576760, 0.6207985133252425, 0.667629, 0.581209),
  ASW_label_batch = c(0.828834, 0.8184154199327823, 0.850226, 0.890412),
  isolated_label_F1 = c(0.135922, 0.12962962962962962, 0.148936, 0.452863),
  graph_conn = c(0.964067, 0.983395565536911, 0.986879, 0.961202)
)
# 保存为CSV
write.csv(df, file = "integration_scores.csv")

library(fmsb)
# 需要最大最小行
radar_df <- rbind(
  max = rep(1, ncol(df)),  # 假设分数上限都是1
  min = rep(0, ncol(df)),  # 分数下限0
  df
)
# 雷达图
radarchart(radar_df, 
           axistype = 1, 
           pcol = c("#F9ED1D", "#9EF8EE", "#B2DBB9", "#F24405"), 
           plwd = 2,
           plty = 1,
           title = "Integration Scores Radar Chart")
legend("topright", legend = rownames(df), col = c("#F9ED1D", "#9EF8EE", "#B2DBB9", "#F24405"), lty = 1, lwd = 2)
