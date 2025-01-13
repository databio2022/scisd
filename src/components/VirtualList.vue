<template>
  <RecycleScroller
    class="scroller"
    :items="list"
    :item-size="24"
    v-slot="{ item }"
    :key="list.length"
  >
    <div class="user" @click="selectGene(item)">
      {{ item }}
    </div>
  </RecycleScroller>
</template>

<script>
import 'vue-virtual-scroller/dist/vue-virtual-scroller.css'
export default {
  props: {
    list: Array
  },
  methods: {
    selectGene(gene) {
      this.$emit('gene-selected', gene)
    }
  }
}
</script>

<style scoped>
.scroller {
  height: 100px;
}

.user {
  height: 24px;
  padding: 0 12px;
  display: flex;
  align-items: center;
}
</style>
