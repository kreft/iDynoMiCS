
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

package utils;


import java.io.File;
//import java.util.ArrayList;
//import java.util.Arrays;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;

import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
//import org.jfree.data.xy.XYDataset;
//import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;


@SuppressWarnings("serial")
public class Chart extends ApplicationFrame {

	private XYSeriesCollection[] _dataSource;
	private String               _resultPath;
	private CombinedDomainXYPlot parentPlot;
	private JFreeChart           chart;
	private ChartPanel           panel;

	public Chart(final String title) {

		super(title);

		parentPlot = new CombinedDomainXYPlot();

		parentPlot.setGap(10.0);
		parentPlot.setOrientation(PlotOrientation.VERTICAL);

		chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, parentPlot, true);
		panel = new ChartPanel(chart, true, true, true, false, true);
		panel.setPreferredSize(new java.awt.Dimension(1024, 768));

		setContentPane(panel);
	}

	public void setPath(String resultPath) {
		_resultPath = resultPath;
	}

	public void init(XYSeriesCollection[] allGraphSets, String[] xLegend, String[] yLegend) {
		_dataSource = allGraphSets;
		XYPlot subplot;
		NumberAxis rangeAxisX, rangeAxisY;
		XYItemRenderer renderer;

		for (int iPlot = 0; iPlot<_dataSource.length; iPlot++) {
			renderer = new StandardXYItemRenderer();
			rangeAxisX = new NumberAxis(xLegend[iPlot]);
			rangeAxisY = new NumberAxis(yLegend[iPlot]);
			subplot = new XYPlot(allGraphSets[iPlot], rangeAxisX, rangeAxisY, renderer);
			subplot.setRangeAxisLocation(AxisLocation.BOTTOM_OR_LEFT);
			parentPlot.add(subplot, 1);
		}
	}

	/**
	 * 
	 * @param iChart
	 * @param iSeries
	 * @param x
	 * @param y
	 */
	public void updateChart(int iChart, int iSeries, double x, double y) {
		_dataSource[iChart].getSeries(iSeries).add(x, y);
	}

	/**
	 * 
	 * @param iChart
	 * @param iSeries
	 * @param x
	 * @param y
	 */
	public void updateChart(int iChart, int iSeries, double[] x, double[] y) {
		_dataSource[iChart].getSeries(iSeries).clear();
		for (int iPoint = 0; iPoint<x.length; iPoint++)
			_dataSource[iChart].getSeries(iSeries).add(x[iPoint], y[iPoint]);

	}

	public void repaintAndSave() {
		// pack();
		this.repaint();
		try {
			ChartUtilities.saveChartAsPNG(new File(_resultPath+"solute.png"), chart, 1024, 768);
		} catch (Exception e) {
		}

		// _frm1.repaint();
		// _frm2.repaint();
		/*
		 * for(Object aPlot:plot.getSubplots()){ ((XYPlot) aPlot). }
		 * 
		 * ChartUtilities.saveChartAsPNG(new
		 * File(_resultPath+"solute.png"),((XYPlot) aPlot) , 1024, 768); } try {
		 * ChartUtilities.saveChartAsPNG(new File(_resultPath+"solute.png"), ,
		 * 1024, 768); ChartUtilities.saveChartAsPNG(new
		 * File(_resultPath+"species.png"), _chartSpecies, 1024, 768); } catch
		 * (Exception e) { }
		 */
	}
}
