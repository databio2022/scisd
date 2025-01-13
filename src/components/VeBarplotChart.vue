<template>
  <div ref="plotContainer"></div>
</template>

<script>
import { ref, onMounted, watch } from 'vue'
import Plotly from 'plotly.js-dist-min'

export default {
  name: 'BarPlotComponent',
  props: {
    data: {
      type: Array,
      required: true
    },
    layout: {
      type: Object,
      default: () => ({})
    }
  },
  setup(props) {
    const plotContainer = ref(null)

    onMounted(() => {
      Plotly.newPlot(plotContainer.value, props.data, props.layout)
    })

    return {
      plotContainer
    }
  }
}
</script>
