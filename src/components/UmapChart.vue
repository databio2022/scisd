<template>
  <div ref="plotlyChart"></div>
</template>

<script>
import { onMounted, ref } from 'vue'
import Plotly from 'plotly.js-dist-min'
//import { UMAP } from 'umap-js'

export default {
  props: {
    data: {
      type: Array,
      required: true
    },
    title: {
      type: String,
      required: true
    }
  },
  setup(props) {
    const plotlyChart = ref(null)

    const drawUmap = (umapData) => {
      const trace = {
        x: umapData.map((item) => item.UMAP_1),
        y: umapData.map((item) => item.UMAP_2),
        mode: 'markers',
        type: 'scatter',
        marker: {
          color: umapData.map((item) => item.x),
          colorscale: 'Reds' // 使用灰色到红色的渐变色
        }
      }

      const layout = {
        title: props.title,
        xaxis: {
          title: 'UMAP_1',
          showgrid: false,
          zeroline: false,
          linewidth: 1,
          linecolor: 'black',
          ticks: 'outside'
        },
        yaxis: {
          title: 'UMAP_2',
          showgrid: false,
          zeroline: false,
          linewidth: 1,
          linecolor: 'black',
          ticks: 'outside'
        },
        autosize: true
      }

      Plotly.newPlot(plotlyChart.value, [trace], layout)
    }

    onMounted(() => {
      drawUmap(props.data)
    })

    return {
      plotlyChart
    }
  }
}
</script>
