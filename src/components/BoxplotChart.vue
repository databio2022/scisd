<template>
  <div ref="plotContainer"></div>
</template>

<script>
import { ref, onMounted } from 'vue'
import Plotly from 'plotly.js-dist-min'

export default {
  name: 'BoxPlotComponent',
  props: {
    data: {
      type: Array,
      required: true
    },
    gene: {
      type: String,
      required: true
    },
    layout: {
      type: Object,
      default: () => ({})
    },
    config: {
      type: Object,
      default: () => ({})
    }
  },
  setup(props) {
    const plotContainer = ref(null)

    onMounted(() => {
      const traces = props.data.map((item) => ({
        y: item.y,
        x: item.x,
        name: item.name,
        type: 'box',
        points: 'none'
      }))

      const layout = {
        title: props.gene + ': the expression profile across cell types',
        yaxis: {
          zeroline: false,
          automargin: true,
        },
        xaxis:{automargin: true,},
        showlegend: false,
        
      }

      Plotly.newPlot(plotContainer.value, traces, layout)
    })

    return {
      plotContainer
    }
  }
}
</script>
